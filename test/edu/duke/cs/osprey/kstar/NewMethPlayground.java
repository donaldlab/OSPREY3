package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.SVG;
import edu.duke.cs.osprey.tools.SVGPlot;
import edu.duke.cs.osprey.tools.TimeTools;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;


public class NewMethPlayground {

	public static final MathContext decimalPrecision = new MathContext(64, RoundingMode.HALF_UP);
	public static final BoltzmannCalculator bcalc = new BoltzmannCalculator(decimalPrecision);

	public static void main(String[] args) {

		File ematFile = new File("emat.dat");
		File confDBFile = new File("newMeth.conf.db");

		Molecule mol = PDBIO.readResource("/1CC8.ss.pdb");

		// TODO: more sequences!
		// make a simple conf space with two sequences
		Strand ligand = new Strand.Builder(mol)
			.setResidues("A2", "A10")
			.build();
		ligand.flexibility.get("A5").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "THR", "PHE", "ARG").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A6").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "THR", "PHE", "ARG").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A7").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "THR", "PHE", "ARG").addWildTypeRotamers().setContinuous();

		Strand target = new Strand.Builder(mol)
			.setResidues("A11", "A72")
			.build();
		target.flexibility.get("A11").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("A12").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("A13").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.addStrand(target)
			.build();

		// calc the emat
		EnergyMatrix emat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()
		) {
			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
			emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(ematFile)
				.build()
				.calcEnergyMatrix();
		}

		final int maxNumBestSequences = 40;
		final int maxNumConfsUpperBounded = 1;

		clearDB(confDBFile);
		calculateUpperBounds(confSpace, emat, confDBFile, maxNumBestSequences);

		calculateLowerBounds(confSpace, emat, confDBFile, maxNumBestSequences, maxNumConfsUpperBounded);
		analyze(confSpace, confDBFile, maxNumBestSequences);
	}

	private static void clearDB(File confDBFile) {
		if (confDBFile.exists()) {
			confDBFile.delete();
		}
	}

	private static void calculateUpperBounds(SimpleConfSpace confSpace, EnergyMatrix emat, File confDBFile, int maxNumBestSequences) {

		new ConfDB(confSpace, confDBFile).use((db) -> {

			ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.setShowProgress(true)
				.build();

			BigInteger numSequences = confSpace.calcNumSequences();
			BigInteger numConfs = astar.getNumConformations();
			log("total A* sequences: %s", formatBig(numSequences));
			log("total A* confs: %s", formatBig(numConfs));

			// calc the max number of confs for any sequence
			BigInteger maxNumConfsForSequence = BigInteger.ZERO;
			for (SimpleConfSpace.Position pos : confSpace.positions) {

				// take the max res confs of any res type
				Map<String,Integer> counts = new HashMap<>();
				for (SimpleConfSpace.ResidueConf resConf : pos.resConfs) {
					Integer count = counts.get(resConf.template.name);
					if (count == null) {
						count = 1;
					} else {
						count++;
					}
					counts.put(resConf.template.name, count);
				}
				BigInteger maxNumResConfs = BigInteger.valueOf(counts.values().stream()
					.max(Comparator.comparingInt((a) -> a))
					.get());

				if (MathTools.isZero(maxNumConfsForSequence)) {
					maxNumConfsForSequence = maxNumResConfs;
				} else {
					maxNumConfsForSequence = maxNumConfsForSequence.multiply(maxNumResConfs);
				}
			}
			log("max confs for any sequence: %s", formatBig(maxNumConfsForSequence));

			log("computing conf energy lower bounds...");

			// keep a list of sequences in order by descending pfLBoUBs, and descending pfUBs
			Map<Sequence,SequenceInfo> infosBySequence = new HashMap<>();
			TreeSet<SequenceInfo> infosByPfuncLBoUB = new TreeSet<>(Comparator.comparing((SequenceInfo info) -> info.pfuncUpperBoundSampled).reversed());
			TreeSet<SequenceInfo> infosByPfuncUB = new TreeSet<>(Comparator.comparing((SequenceInfo info) -> info.pfuncUpperBound).reversed());

			// compute the top K sequences by pfUBs
			List<SequenceInfo> bestInfos = new ArrayList<>();
			while (true) {

				// get the next conf and energy lower bound
				ConfSearch.ScoredConf conf = astar.nextConf();
				if (conf == null) {
					log("%d sequences found, but no other sequences exist, so we're done", infosBySequence.size());
					break;
				}

				// check the boltzmann-weighted energy
				BigDecimal weightedLowerEnergy = bcalc.calc(conf.getScore());
				if (MathTools.isZero(weightedLowerEnergy)) {
					log("conf space ran out of low energies. There's no more useful information to be learned, so we have to stop");
					break;
				}

				// update the conf db
				ConfDB.SequenceDB sdb = db.getSequence(confSpace.makeSequenceFromConf(conf));
				sdb.setLowerBound(
					conf.getAssignments(),
					conf.getScore(),
					TimeTools.getTimestampNs()
				);
				sdb.updateLowerEnergyOfUnsampledConfs(conf.getScore());

				// get (or create) the sequence info
				SequenceInfo info = infosBySequence.get(sdb.sequence);
				if (info == null) {
					info = new SequenceInfo(sdb.sequence);
					infosBySequence.put(sdb.sequence, info);

					log("discovered %d/%d unique sequences", infosBySequence.size(), numSequences);
				}

				// update pfunc upper bound for this sequence
				info.updatePfuncUpperBound(conf.getScore());

				// re-sort sequences
				for (TreeSet<SequenceInfo> infos : Arrays.asList(infosByPfuncLBoUB, infosByPfuncUB)) {
					infos.remove(info);
					infos.add(info);
				}

				// compute an upper bound on the pfunc for all unsampled sequences
				BigDecimal otherSequencesPfuncUB = weightedLowerEnergy.multiply(new BigDecimal(maxNumConfsForSequence), decimalPrecision);

				/* TEMP
				log("%4d/%d   sequence: %s   ceLB: %.3f   wceLB: %s   pfUB: [%s,%s]   pfUB for unsampled sequences: %s",
					++numConfsLowerBounded, numConfs,
					sdb.sequence.toString(),
					conf.getScore(), formatBig(weightedLowerEnergy),
					formatBig(info.pfuncUpperBoundSampled), formatBig(info.pfuncUpperBound),
					formatBig(otherSequencesPfuncUB)
				);
				*/

				// grab the best K sequences by pfLBoUB
				bestInfos.clear();
				for (SequenceInfo otherInfo : infosByPfuncLBoUB) {
					bestInfos.add(otherInfo);
					if (bestInfos.size() >= maxNumBestSequences) {
						break;
					}
				}

				Supplier<Boolean> foundBestSequences = () -> {

					// could an unsampled sequence have a higher pfUB than a best pfLBoUB?
					for (SequenceInfo bestInfo : bestInfos) {
						if (MathTools.isGreaterThanOrEqual(otherSequencesPfuncUB, bestInfo.pfuncUpperBoundSampled)) {
							return false;
						}
					}

					// does the best remaining sequence by pfUB have a higher pfUB?
					SequenceInfo bestRemaining = infosByPfuncUB.stream()
						.filter((otherInfo) -> !bestInfos.contains(otherInfo))
						.findFirst()
						.orElse(null);
					if (bestRemaining == null) {
						return false;
					}
					for (SequenceInfo bestInfo : bestInfos) {
						if (MathTools.isGreaterThanOrEqual(bestRemaining.pfuncUpperBound, bestInfo.pfuncUpperBoundSampled)) {
							return false;
						}
					}

					// the best K really are the best! =D
					return true;
				};

				// are we done yet?
				if (foundBestSequences.get()) {
					break;
				}
			}

			// show the results
			for (SequenceInfo info : bestInfos) {
				System.out.println(info.toString());
			}

			// plot the bounds on pfunc upper bounds
			SVG svg = new SVG();

			SVG.StyleClass wildtypeIntervalStyle = svg.makeStyleClass("wildtype-interval-style");
			wildtypeIntervalStyle.setStrokeColor(0x66cc55);
			wildtypeIntervalStyle.setStrokeWidth(0.5);

			SVGPlot.Intervals intervals = new SVGPlot.Intervals();
			intervals.intervalWidth = 2.0;
			for (SequenceInfo info : bestInfos) {
				SVGPlot.Intervals.Interval interval = intervals.addInterval(
					Math.log10(info.pfuncUpperBoundSampled.doubleValue()),
					Math.log10(info.pfuncUpperBound.doubleValue())
				);
				interval.id = info.sequence.toString();
				if (info.sequence.isWildType()) {
					interval.extraStyle = wildtypeIntervalStyle;
				}
			}
			intervals.draw(svg);
			intervals.setBounds(svg, 10, 16);
			svg.finish().write(new File("pfuncUpperBounds.svg"));

		}); // db
	}

	private static void calculateLowerBounds(SimpleConfSpace confSpace, EnergyMatrix emat, File confDBFile, int maxNumBestSequences, int maxNumConfsUpperBounded) {

		new ConfDB(confSpace, confDBFile).use((db) -> {

			// get all our sequence info, sorted by descending pfUB
			TreeSet<SequenceInfo> infosByPfuncUB = new TreeSet<>(Comparator.comparing((SequenceInfo info) -> info.pfuncUpperBound).reversed());
			for (Sequence sequence : db.getSequences()) {
				SequenceInfo info = new SequenceInfo(sequence);
				info.readPfuncUpperBoundFromDB(db.getSequence(sequence));
				infosByPfuncUB.add(info);
			}

			new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
				.setParallelism(Parallelism.makeCpu(4))
				.use((ecalc) -> {

					ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

					// estimate pfLB for the top K sequences by pfUB
					int numSequencesLowerBounded = 0;
					for (SequenceInfo info : infosByPfuncUB) {

						log("estimating pfunc lower bound for sequence %s", info.sequence);

						ConfDB.SequenceDB sdb = db.getSequence(info.sequence);
						info.calcPfuncLowerBound(confEcalc, sdb, maxNumConfsUpperBounded);

						if (++numSequencesLowerBounded >= maxNumBestSequences) {
							break;
						}
					}

				}); // ecalc
		});
	}

	private static void analyze(SimpleConfSpace confSpace, File confDBFile, int maxNumBestSequences) {

		new ConfDB(confSpace, confDBFile).use((db) -> {

			// get all our sequence info, sorted by descending pfUB
			TreeSet<SequenceInfo> infosByPfuncUB = new TreeSet<>(Comparator.comparing((SequenceInfo info) -> info.pfuncUpperBound).reversed());
			for (Sequence sequence : db.getSequences()) {
				SequenceInfo info = new SequenceInfo(sequence);
				ConfDB.SequenceDB sdb = db.getSequence(sequence);
				info.readPfuncUpperBoundFromDB(sdb);
				info.readPfuncLowerBoundFromDB(sdb);
				infosByPfuncUB.add(info);
			}

			// calculate the ranking indices, RL and RU
			// meaning, a sequence is in the top K sequences for K = RL(pfLB)
			// and a sequence is NOT in the top K sequences for K = RU(pfUB)
			TreeMap<BigDecimal,Integer> rankIndexLower = new TreeMap<>();
			TreeMap<BigDecimal,Integer> rankIndexUpper = new TreeMap<>();
			for (SequenceInfo info : infosByPfuncUB) {
				rankIndexLower.put(info.pfuncLowerBound, rankIndexLower.size() + 2);
				rankIndexUpper.put(info.pfuncUpperBound, rankIndexUpper.size());
			}
			log("indexed sequences to top %d", infosByPfuncUB.size());

			Function<BigDecimal,Integer> getRankUpper = (pfuncLB) -> {
				Map.Entry<BigDecimal,Integer> entry = rankIndexUpper.floorEntry(pfuncLB);
				if (entry != null) {
					return entry.getValue();
				}
				return null;
			};

			Function<BigDecimal,Integer> getRankLower = (pfuncUB) -> {
				Map.Entry<BigDecimal,Integer> entry = rankIndexLower.ceilingEntry(pfuncUB);
				if (entry != null) {
					return entry.getValue();
				}
				return 1;
			};

			// show the ranked sequences
			for (SequenceInfo info : infosByPfuncUB) {
				int rankLower = getRankLower.apply(info.pfuncUpperBound);
				Integer rankUpper = getRankUpper.apply(info.pfuncLowerBound);
				if (rankUpper != null) {
					log("%s is ranked in [%d,%d]", info.sequence, rankLower, rankUpper);
				}
			}

			log("pfunc cutoff for top %d sequences is %s", infosByPfuncUB.size() - 1, formatBig(infosByPfuncUB.last().pfuncUpperBound));

			// plot the pfunc bounds for the best K sequences
			SVG svg = new SVG();

			SVG.StyleClass wildtypeIntervalStyle = svg.makeStyleClass("wildtype-interval-style");
			wildtypeIntervalStyle.setStrokeColor(0x66cc55);
			wildtypeIntervalStyle.setStrokeWidth(0.5);

			SVGPlot.Intervals intervals = new SVGPlot.Intervals();
			intervals.intervalWidth = 2.0;
			int numSequences = 0;
			for (SequenceInfo info : infosByPfuncUB) {
				SVGPlot.Intervals.Interval interval = intervals.addInterval(
					Math.max(0, Math.log10(info.pfuncLowerBound.doubleValue())),
					Math.log10(info.pfuncUpperBound.doubleValue())
				);
				interval.id = String.format("%s: [%s,%s]",
					info.sequence.toString(),
					formatBig(info.pfuncLowerBound),
					formatBig(info.pfuncUpperBoundSampled)
				);
				if (info.sequence.isWildType()) {
					interval.extraStyle = wildtypeIntervalStyle;
				}
				if (++numSequences >= maxNumBestSequences) {
					break;
				}
			}
			intervals.draw(svg);
			intervals.setBounds(svg, 10, 16);

			// show the lowest pfUB of any sequence we know
			{
				SequenceInfo info = infosByPfuncUB.last();
				double y = Math.log10(info.pfuncUpperBound.doubleValue());

				SVG.StyleClass dottedLineStyle = svg.makeStyleClass("dottedLineStyle");
				dottedLineStyle.setStrokeColor(0xcccccc);
				dottedLineStyle.setStrokeWidth(0.2);
				dottedLineStyle.setStrokeDashArray(1, 1);
				svg.makeLine(
						intervals.xmin + intervals.intervalSpacing, y,
						intervals.xmax, y
					)
					.setStyleClasses(dottedLineStyle)
					.setId(String.format("cutoff for top %d sequences", infosByPfuncUB.size() - 1))
					.draw();
			}

			svg.finish().write(new File("pfuncs.svg"));
		});
	}

	public static class SequenceInfo {

		public final Sequence sequence;

		public BigInteger numConfs;
		public long numConfsLowerBounded;
		public long numConfsUpperBounded;

		public BigDecimal pfuncUpperBoundSampled = null;
		public BigDecimal pfuncUpperBoundUnsampled = null;
		public BigDecimal pfuncUpperBound = null;
		public BigDecimal pfuncLowerBound = null;

		public SequenceInfo(Sequence sequence) {

			this.sequence = sequence;

			// count the number of confs
			numConfs = BigInteger.ZERO;
			for (SimpleConfSpace.Position pos : sequence.confSpace.positions) {
				String resType = sequence.get(pos);
				BigInteger numResConfs = BigInteger.valueOf(pos.resConfs.stream()
					.filter((resConf) -> resConf.template.name.equals(resType))
					.count());

				if (MathTools.isZero(numConfs)) {
					numConfs = numResConfs;
				} else {
					numConfs = numConfs.multiply(numResConfs);
				}
			}

			numConfsLowerBounded = 0;
			numConfsUpperBounded = 0;
		}

		public void readPfuncUpperBoundFromDB(ConfDB.SequenceDB sdb) {

			pfuncUpperBoundSampled = BigDecimal.ZERO;
			numConfsLowerBounded = 0;

			for (ConfDB.Conf conf : sdb) {
				if (conf.lower != null) {
					pfuncUpperBoundSampled = MathTools.bigAdd(
						pfuncUpperBoundSampled,
						bcalc.calc(conf.lower.energy),
						decimalPrecision
					);
					numConfsLowerBounded++;
				}
			}

			pfuncUpperBoundUnsampled = bcalc.calc(sdb.getLowerEnergyOfUnsampledConfs()).multiply(new BigDecimal(getNumConfsNotLowerBounded()), decimalPrecision);
			pfuncUpperBound = pfuncUpperBoundSampled.add(pfuncUpperBoundUnsampled, decimalPrecision);
		}

		public void readPfuncLowerBoundFromDB(ConfDB.SequenceDB sdb) {

			numConfsUpperBounded = 0;
			pfuncLowerBound = BigDecimal.ZERO;

			for (ConfDB.Conf conf : sdb) {
				if (conf.upper != null) {
					pfuncLowerBound = MathTools.bigAdd(
						pfuncLowerBound,
						bcalc.calc(conf.upper.energy),
						decimalPrecision
					);
					numConfsUpperBounded++;
				}
			}
		}

		/**
		 * @param lowerEnergy must always increase between successive calls
		 */
		public void updatePfuncUpperBound(double lowerEnergy) {

			// init state if this the first time update has been called
			if (pfuncUpperBound == null) {
				pfuncUpperBoundSampled = BigDecimal.ZERO;
				numConfsLowerBounded = 0;
				pfuncUpperBoundUnsampled = MathTools.BigPositiveInfinity;
				pfuncUpperBound = MathTools.BigPositiveInfinity;
			}

			BigDecimal weightedLowerEnergy = bcalc.calc(lowerEnergy);

			pfuncUpperBoundSampled = MathTools.bigAdd(
				pfuncUpperBoundSampled,
				weightedLowerEnergy,
				decimalPrecision
			);
			numConfsLowerBounded++;

			pfuncUpperBoundUnsampled = weightedLowerEnergy.multiply(new BigDecimal(getNumConfsNotLowerBounded()), decimalPrecision);
			pfuncUpperBound = pfuncUpperBoundSampled.add(pfuncUpperBoundUnsampled, decimalPrecision);
		}

		public void calcPfuncLowerBound(ConfEnergyCalculator confEcalc, ConfDB.SequenceDB sdb, int maxNumConfs) {

			// sort confs by lower bound
			List<ConfDB.Conf> confs = new ArrayList<>();
			for (ConfDB.Conf conf : sdb) {
				confs.add(conf);
			}
			confs.sort(Comparator.comparingDouble((conf) -> conf.lower.energy));

			// TODO: allow time-based limits instead of conf-based limits?
			// TODO: use smart early exit of some kind?

			// TODO: skip energies we've already calculated

			// estimate the pfunc lower bound
			pfuncLowerBound = BigDecimal.ZERO;
			maxNumConfs = Math.min(maxNumConfs, confs.size());
			for (ConfDB.Conf conf : confs.subList(0, maxNumConfs)) {
				RCTuple frag = new RCTuple(conf.assignments);
				confEcalc.calcEnergyAsync(frag, (epmol) -> {

					// NOTE: this is on a listener thread, separate from the main thread

					sdb.setUpperBound(
						conf.assignments,
						epmol.energy,
						TimeTools.getTimestampNs()
					);

					pfuncLowerBound = MathTools.bigAdd(
						pfuncLowerBound,
						bcalc.calc(epmol.energy),
						decimalPrecision
					);

					numConfsUpperBounded++;
				});
			}

			confEcalc.ecalc.tasks.waitForFinish();
		}

		public BigInteger getNumConfsNotLowerBounded() {
			return numConfs.subtract(BigInteger.valueOf(numConfsLowerBounded));
		}

		@Override
		public String toString() {
			StringBuilder buf = new StringBuilder();
			log(buf, "sequence: %s", sequence);
			log(buf, "\ttotal confs:                   %s", formatBig(numConfs));
			log(buf, "\tnum confs lower bounded:       %d", numConfsLowerBounded);
			log(buf, "\tnum confs upper bounded:       %d", numConfsUpperBounded);
			if (pfuncUpperBound != null) {
				log(buf, "\tpfunc upper bound:             %s", formatBig(pfuncUpperBound));
				log(buf, "\t\tfrom sampled:                    %s", formatBig(pfuncUpperBoundSampled));
				log(buf, "\t\tfrom unsampled:                  %s", formatBig(pfuncUpperBoundUnsampled));
				log(buf, "\t\tpercent sampled:                 %.4f", MathTools.bigDivide(pfuncUpperBoundSampled, pfuncUpperBound, decimalPrecision).doubleValue()*100.0f);
			}
			if (pfuncLowerBound != null) {
				log(buf, "\tpfunc lower bound:             %s", formatBig(pfuncLowerBound));
				double logUncertainty = Math.log10(pfuncUpperBound.doubleValue()) - Math.log10(pfuncLowerBound.doubleValue());
				log(buf, "\tpfunc log uncertainty:         %.2f", logUncertainty);
			}
			return buf.toString();
		}
	}

	private static void log(String format, Object ... args) {
		System.out.println(String.format(format, args));
	}

	private static void log(StringBuilder buf, String format, Object ... args) {
		buf.append(String.format(format + "\n", args));
	}

	private static String formatBig(BigInteger i) {
		if (i.compareTo(BigInteger.valueOf(1000000)) < 0) {
			return String.format("%s", i);
		} else {
			return String.format("%e", i.doubleValue());
		}
	}

	private static String formatBig(BigDecimal f) {
		return String.format("%e (%.2f)", f.doubleValue(), Math.log10(f.doubleValue()));
	}
}
