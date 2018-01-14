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

		// TODO: more sequences!
		// make a simple conf space with two sequences
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "THR", "PHE", "ARG").addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A6").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "THR", "PHE", "ARG").addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A7").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "THR", "PHE", "ARG").addWildTypeRotamers().setContinuous();

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
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

		final long maxNumConfsLowerBounded = 10000L;
		final int maxNumBestSequences = 5;
		final int maxNumConfsUpperBounded = 1000;

		clearDB(confDBFile);
		calculateUpperBounds(confSpace, emat, confDBFile, maxNumBestSequences);
		//calculateFixed(confSpace, emat, confDBFile, maxNumConfsLowerBounded, maxNumBestSequences, maxNumConfsUpperBounded);
		//analyze(confSpace, confDBFile, maxNumBestSequences);
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
			long numConfsLowerBounded = 0L;
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

					log("found %d/%d sequences", infosBySequence.size(), numSequences);
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

				// TEMP
				log("%4d/%d   sequence: %s   ceLB: %.3f   wceLB: %s   pfUB: [%s,%s]   pfUB for unsampled sequences: %s",
					++numConfsLowerBounded, numConfs,
					sdb.sequence.toString(),
					conf.getScore(), formatBig(weightedLowerEnergy),
					formatBig(info.pfuncUpperBoundSampled), formatBig(info.pfuncUpperBound),
					formatBig(otherSequencesPfuncUB)
				);

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
			SVGPlot.Intervals intervals = new SVGPlot.Intervals();
			for (SequenceInfo info : bestInfos) {
				intervals.addInterval(
					Math.log10(info.pfuncUpperBoundSampled.doubleValue()),
					Math.log10(info.pfuncUpperBound.doubleValue())
				);
			}
			intervals.draw(svg);
			intervals.setBounds(svg, 10, 16);
			svg.finish().write(new File("pfuncUpperBounds.svg"));

		}); // db
	}

	private static void calculateFixed(SimpleConfSpace confSpace, EnergyMatrix emat, File confDBFile, long maxNumConfsLowerBounded, int maxNumBestSequences, int maxNumConfsUpperBounded) {

		new ConfDB(confSpace, confDBFile).use((db) -> {

			new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
				.setParallelism(Parallelism.makeCpu(4))
				.use((ecalc) -> {

					ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

					ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
						.setTraditional()
						.setShowProgress(true)
						.build();

					BigInteger numSequences = confSpace.calcNumSequences();
					log("total A* sequences: %s", formatBig(numSequences));
					log("total A* confs: %s", formatBig(astar.getNumConformations()));

					// compute lower bounds for the first K conformations in the A* tree
					long numConfsLowerBounded = 0L;
					while (true) {

						ConfSearch.ScoredConf conf = astar.nextConf();
						if (conf == null) {
							break;
						}

						ConfDB.SequenceDB sdb = db.getSequence(confSpace.makeSequenceFromConf(conf));

						sdb.setLowerBound(
							conf.getAssignments(),
							conf.getScore(),
							TimeTools.getTimestampNs()
						);
						sdb.updateLowerEnergyOfUnsampledConfs(conf.getScore());

						if (++numConfsLowerBounded >= maxNumConfsLowerBounded) {
							break;
						}
					}

					log("Finished lower bounding conf energies");

					// sort sequences by pfunc descending upper bounds
					List<SequenceInfo> sequencesByPfuncUpperBound = new ArrayList<>();
					for (Sequence sequence : db.getSequences()) {
						SequenceInfo info = new SequenceInfo(sequence);
						info.readPfuncUpperBoundFromDB(db.getSequence(sequence));
						sequencesByPfuncUpperBound.add(info);
					}
					sequencesByPfuncUpperBound.sort((a, b) -> b.pfuncUpperBound.compareTo(a.pfuncUpperBound));

					// pick the best K sequences by pfunc upper bounds
					int numBestSequences = Math.min(maxNumBestSequences, sequencesByPfuncUpperBound.size());
					List<SequenceInfo> bestSequences = sequencesByPfuncUpperBound.subList(0, numBestSequences);
					log("%d/%s best sequences by PFunc upper bounds:", maxNumBestSequences, formatBig(numSequences));
					for (SequenceInfo info : bestSequences) {
						log("Best Sequences: %s -> %e", info.sequence, info.pfuncUpperBound);
					}

					// estimate lower bounds for the best sequences using at most K confs
					for (SequenceInfo info : bestSequences) {
						log("Estimating PFunc lower bound for %s ...", info.sequence);
						info.calcPfuncLowerBound(confEcalc, db.getSequence(info.sequence), maxNumConfsUpperBounded);
					}

				}); // ecalc
		});
	}

	private static void analyze(SimpleConfSpace confSpace, File confDBFile, int maxNumBestSequences) {

		new ConfDB(confSpace, confDBFile).use((db) -> {

			// sort sequences by pfunc descending upper bounds
			List<SequenceInfo> sequencesByPfuncUpperBound = new ArrayList<>();
			for (Sequence sequence : db.getSequences()) {
				SequenceInfo info = new SequenceInfo(sequence);
				info.readPfuncUpperBoundFromDB(db.getSequence(sequence));
				sequencesByPfuncUpperBound.add(info);
			}
			sequencesByPfuncUpperBound.sort((a, b) -> b.pfuncUpperBound.compareTo(a.pfuncUpperBound));
			int numBestSequences = Math.min(maxNumBestSequences, sequencesByPfuncUpperBound.size());
			List<SequenceInfo> bestSequences = sequencesByPfuncUpperBound.subList(0, numBestSequences);

			// calc the pfunc lower bounds too
			for (SequenceInfo info : bestSequences) {
				info.readPfuncLowerBoundFromDB(db.getSequence(info.sequence));
			}

			// plot the pfunc bounds
			SVG svg = new SVG();
			SVGPlot.Intervals intervals = new SVGPlot.Intervals();
			for (SequenceInfo info : bestSequences) {
				intervals.addInterval(
					Math.log10(info.pfuncLowerBound.doubleValue()),
					Math.log10(info.pfuncUpperBound.doubleValue())
				);
			}
			intervals.draw(svg);
			intervals.setBounds(svg, 10, 10);
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

			// TODO: use time-based limits instead of conf-based limits?
			// TODO: use smart early exit of some kind?

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
