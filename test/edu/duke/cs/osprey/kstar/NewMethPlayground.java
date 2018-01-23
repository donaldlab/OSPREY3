package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.pruning.AStarSequencePruner;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.*;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;


public class NewMethPlayground {

	public static final MathContext decimalPrecision = new MathContext(64, RoundingMode.HALF_UP);
	public static final BoltzmannCalculator bcalc = new BoltzmannCalculator(decimalPrecision);

	public static void main(String[] args) {

		File complexEmatFile = new File("emat.complex.dat");
		File ligandEmatFile = new File("emat.ligand.dat");
		File complexConfDBFile = new File("conf.complex.db");
		File ligandConfDBFile = new File("conf.ligand.db");

		// try JJ's design
		Molecule mol = PDBIO.readResource("/gp120/gp120SRVRC26.09SR.pdb");

		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder()
			.addMoleculeForWildTypeRotamers(mol)
			.addTemplates(FileTools.readResource("/gp120/all_nuc94_and_gr.in"))
			.addTemplateCoords(FileTools.readResource("/gp120/all_amino_coords.in"))
			.addRotamers(FileTools.readResource("/gp120/GenericRotamers.dat"))
			.build();

		Strand ligand = new Strand.Builder(mol)
			.setResidues("H1792", "L2250")
			.setTemplateLibrary(templateLib)
			.build();
		ligand.flexibility.get("H1901").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1904").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIS", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1905").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIS", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1906").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1907").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1908").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		Strand target = new Strand.Builder(mol)
			.setResidues("F379", "J1791")
			.setTemplateLibrary(templateLib)
			.build();
		target.flexibility.get("G973").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G977").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G978").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G979").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G980").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("J1448").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		SimpleConfSpace complexConfSpace = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.addStrand(target)
			.build();

		SimpleConfSpace ligandConfSpace = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();

		SimpleConfSpace targetConfSpace = new SimpleConfSpace.Builder()
			.addStrand(target)
			.build();

		// calc the emats
		EnergyMatrix complexEmat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(complexConfSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()
		) {
			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(complexConfSpace, ecalc).build();
			complexEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(complexEmatFile)
				.build()
				.calcEnergyMatrix();
		}

		// calc the emats
		EnergyMatrix ligandEmat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(ligandConfSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()
		) {
			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(ligandConfSpace, ecalc).build();
			ligandEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(ligandEmatFile)
				.build()
				.calcEnergyMatrix();
		}

		final int maxNumBestSequences = 100;
		final int maxNumConfsPerSequence = 100;
		final int maxNumConfsUpperBounded = 1;
		final double pfuncUBFractionSampled = 0.99;

		// analyze complex
		//clearDB(complexConfDBFile);
		//calculateComplexUpperBounds(complexConfSpace, complexEmat, complexConfDBFile, maxNumBestSequences, maxNumConfsPerSequence);
		//calculateComplexLowerBounds(complexConfSpace, complexEmat, complexConfDBFile, maxNumBestSequences, maxNumConfsUpperBounded);
		//List<SequenceInfo> bestComplexes = analyzeComplexes(complexConfSpace, complexConfDBFile, maxNumBestSequences);

		// analyze ligand for the best complex sequences
		//clearDB(ligandConfDBFile);
		//calculateLigandBounds(ligandConfSpace, ligandEmat, ligandConfDBFile, bestComplexes, pfuncUBFractionSampled, maxNumConfsUpperBounded);
		//analyze(complexConfSpace, complexConfDBFile, maxNumBestSequences, ligandConfSpace, ligandConfDBFile);
	}

	private static void clearDB(File confDBFile) {
		if (confDBFile.exists()) {
			confDBFile.delete();
		}
	}

	private static void calculateComplexUpperBounds(SimpleConfSpace confSpace, EnergyMatrix emat, File confDBFile, int maxNumBestSequences, double maxNumConfsPerSequence) {

		new ConfDB(confSpace, confDBFile).use((db) -> {

			AStarSequencePruner pruner = new AStarSequencePruner(confSpace);

			Supplier<ConfAStarTree> astarFactory = () -> {
				ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
					//.setTraditional()
					.setMPLP(new ConfAStarTree.MPLPBuilder()
						.setUpdater(new EdgeUpdater())
						.setNumIterations(40)
					)
					.setPruner(pruner)
					.setShowProgress(true)
					.build();
				astar.setParallelism(Parallelism.makeCpu(6));
				return astar;
			};

			ConfAStarTree astar = astarFactory.get();

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

				// if this sequence's pfUB is sufficiently estimated, prune all confs for that sequence from the tree
				if (info.numConfsLowerBounded >= maxNumConfsPerSequence) {
					// TEMP
					log("estimated pfUB to %.12f for %s", info.getPfuncUpperBoundFractionSampled(), info.sequence);
					pruner.add(info.sequence);

					// recycle the A* tree
					astar = astarFactory.get();
				}

				// re-sort sequences
				for (TreeSet<SequenceInfo> infos : Arrays.asList(infosByPfuncLBoUB, infosByPfuncUB)) {
					infos.remove(info);
					infos.add(info);
				}

				// compute an upper bound on the pfunc for all unsampled sequences
				BigDecimal otherSequencesPfuncUB = weightedLowerEnergy.multiply(new BigDecimal(maxNumConfsForSequence), decimalPrecision);

				/* TEMP
				log("sequence: %s   ceLB: %.3f   wceLB: %s   pfUB: [%s,%s] %.6f   pfUB for unsampled sequences: %s",
					info.sequence,
					conf.getScore(), formatBig(weightedLowerEnergy),
					formatBig(info.pfuncUpperBoundSampled), formatBig(info.pfuncUpperBound),
					info.getPfuncUpperBoundFractionSampled(),
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
			wildtypeIntervalStyle.priority = 10;
			wildtypeIntervalStyle.setStrokeColor(0x66cc55);
			wildtypeIntervalStyle.setStrokeWidth(0.5);

			SVGPlot.Intervals intervals = new SVGPlot.Intervals();
			intervals.intervalWidth = 2.0;
			for (SequenceInfo info : bestInfos) {
				SVGPlot.Intervals.Interval interval = intervals.addInterval(
					MathTools.log10p1(info.pfuncUpperBoundSampled),
					MathTools.log10p1(info.pfuncUpperBound)
				);
				interval.id = info.sequence.toString();
				if (info.sequence.isWildType()) {
					interval.extraStyle = wildtypeIntervalStyle;
				}
			}
			intervals.axis = intervals.makeAxis();
			intervals.axis.tickFormat = "%.0f";
			intervals.draw(svg);
			intervals.setBounds(svg, 10, 16);
			svg.finish().write(new File("pfunc.complex.upperBounds.svg"));

		}); // db
	}

	private static void calculateComplexLowerBounds(SimpleConfSpace confSpace, EnergyMatrix emat, File confDBFile, int maxNumBestSequences, int maxNumConfsUpperBounded) {

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

	private static List<SequenceInfo> analyzeComplexes(SimpleConfSpace confSpace, File confDBFile, int maxNumBestSequences) {

		List<SequenceInfo> bestSequences = new ArrayList<>();

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

			// get the best K sequences
			int numSequences = 0;
			for (SequenceInfo info : infosByPfuncUB) {
				bestSequences.add(info);
				if (++numSequences >= maxNumBestSequences) {
					break;
				}
			}

			// plot the pfunc bounds for the best K sequences
			SVG svg = new SVG();

			SVG.StyleClass wildtypeIntervalStyle = svg.makeStyleClass("wildtype-interval-style");
			wildtypeIntervalStyle.priority = 10;
			wildtypeIntervalStyle.setStrokeColor(0x66cc55);
			wildtypeIntervalStyle.setStrokeWidth(0.5);

			SVGPlot.Intervals intervals = new SVGPlot.Intervals();
			intervals.intervalWidth = 2.0;
			for (SequenceInfo info : bestSequences) {
				SVGPlot.Intervals.Interval interval = intervals.addInterval(
					MathTools.log10p1(info.pfuncLowerBound),
					MathTools.log10p1(info.pfuncUpperBound)
				);
				interval.id = String.format("%s: [%s,%s]",
					info.sequence.toString(),
					formatBig(info.pfuncLowerBound),
					formatBig(info.pfuncUpperBoundSampled)
				);
				if (info.sequence.isWildType()) {
					interval.extraStyle = wildtypeIntervalStyle;
				}
			}
			intervals.axis = intervals.makeAxis();
			intervals.axis.tickFormat = "%.0f";
			intervals.draw(svg);
			intervals.setBounds(svg, 10, 16);

			// show the lowest pfUB of any sequence we know
			{
				SequenceInfo info = infosByPfuncUB.last();
				double y = MathTools.log10p1(info.pfuncUpperBound);

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

			svg.finish().write(new File("pfuncs.complex.svg"));
		});

		return bestSequences;
	}

	private static Sequence complexToLigandSequence(SimpleConfSpace ligandConfSpace, Sequence complexSequence) {
		Sequence ligandSequence = ligandConfSpace.makeUnassignedSequence();
		for (SimpleConfSpace.Position complexPos : complexSequence.confSpace.positions) {
			SimpleConfSpace.Position ligandPos = ligandConfSpace.getPositionOrNull(complexPos.resNum);
			if (ligandPos != null) {
				ligandSequence.set(ligandPos, complexSequence.get(complexPos));
			}
		}
		assert (ligandSequence.isFullyAssigned());
		return ligandSequence;
	}

	private static void calculateLigandBounds(SimpleConfSpace confSpace, EnergyMatrix emat, File confDBFile, List<SequenceInfo> bestComplexes, double pfuncUBFactionSampled, int maxNumConfsUpperBounded) {

		new ConfDB(confSpace, confDBFile).use((db) -> {

			new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
				.setParallelism(Parallelism.makeCpu(4))
				.use((ecalc) -> {

					ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

					for (SequenceInfo complex : bestComplexes) {

						// convert complex sequence to ligand sequence
						SequenceInfo ligand = new SequenceInfo(complexToLigandSequence(confSpace, complex.sequence));

						ConfDB.SequenceDB sdb = db.getSequence(ligand.sequence);

						log("calculating ligand pfunc for %s ...", complex.sequence);

						ConfSearch.ScoredConf minBoundConf = null;

						// calculate pfunc upper bound
						ConfAStarTree astar = new ConfAStarTree.Builder(emat, ligand.sequence.makeRCs())
							.setTraditional()
							.build();
						while (true) {

							ConfSearch.ScoredConf conf = astar.nextConf();
							if (conf == null) {
								break;
							}

							if (minBoundConf == null) {
								minBoundConf = conf;
							}

							ligand.updatePfuncUpperBound(conf.getScore());

							// update the conf db
							sdb.setLowerBound(
								conf.getAssignments(),
								conf.getScore(),
								TimeTools.getTimestampNs()
							);
							sdb.updateLowerEnergyOfUnsampledConfs(conf.getScore());

							// are we done yet?
							if (ligand.getPfuncUpperBoundFractionSampled() > pfuncUBFactionSampled) {
								break;
							}
						}

						// compute a crude pfunc lower bound
						ligand.calcPfuncLowerBound(confEcalc, sdb, maxNumConfsUpperBounded);

						log("ligand: %s", ligand);
					}
				}); // ecalc
		}); // db
	}

	static class SequenceInfoPair {

		public final SequenceInfo complex;
		public final SequenceInfo ligand;

		public SequenceInfoPair(SequenceInfo complex, SequenceInfo ligand) {
			this.complex = complex;
			this.ligand = ligand;
		}
	}

	private static void analyze(SimpleConfSpace complexConfSpace, File complexConfDBFile, int maxNumBestSequences, SimpleConfSpace ligandConfSpace, File ligandConfDBFile) {

		List<SequenceInfoPair> pairs = new ArrayList<>();

		new ConfDB(complexConfSpace, complexConfDBFile).use((complexDB) -> {

			// get all our sequence info, sorted by descending pfUB
			TreeSet<SequenceInfo> infosByPfuncUB = new TreeSet<>(Comparator.comparing((SequenceInfo info) -> info.pfuncUpperBound).reversed());
			for (Sequence sequence : complexDB.getSequences()) {
				SequenceInfo info = new SequenceInfo(sequence);
				ConfDB.SequenceDB sdb = complexDB.getSequence(sequence);
				info.readPfuncUpperBoundFromDB(sdb);
				info.readPfuncLowerBoundFromDB(sdb);
				infosByPfuncUB.add(info);
			}

			// get the best K complexes by pfUB
			List<SequenceInfo> bestComplexes = new ArrayList<>();
			int numSequences = 0;
			for (SequenceInfo info : infosByPfuncUB) {
				bestComplexes.add(info);
				if (++numSequences >= maxNumBestSequences) {
					break;
				}
			}

			new ConfDB(ligandConfSpace, ligandConfDBFile).use((ligandDB) -> {

				// get the ligand info for those best K sequences
				for (SequenceInfo complex : bestComplexes) {

					SequenceInfo ligand = new SequenceInfo(complexToLigandSequence(ligandConfSpace, complex.sequence));
					ConfDB.SequenceDB sdb = ligandDB.getSequence(ligand.sequence);
					ligand.readPfuncLowerBoundFromDB(sdb);
					ligand.readPfuncUpperBoundFromDB(sdb);

					pairs.add(new SequenceInfoPair(complex, ligand));
				}

			}); // ligand db

		}); // complex db

		// plot the pfunc bounds
		{
			SVG svg = new SVG();

			SVG.StyleClass wildtypeBoxStyle = svg.makeStyleClass("wildtype-box-style");
			wildtypeBoxStyle.priority = 10;
			wildtypeBoxStyle.setStrokeColor(0x66cc55);
			wildtypeBoxStyle.setStrokeWidth(0.5);

			SVGPlot.Boxes boxes = new SVGPlot.Boxes();
			boxes.boxStyle.setNoFill();
			for (SequenceInfoPair pair : pairs) {
				SVGPlot.Boxes.Box box = boxes.addBox(
					MathTools.log10p1(pair.ligand.pfuncLowerBound),
					MathTools.log10p1(pair.ligand.pfuncUpperBound),
					MathTools.log10p1(pair.complex.pfuncLowerBound),
					MathTools.log10p1(pair.complex.pfuncUpperBound)
				);
				box.id = String.format("%s: complex:[%s,%s] ligand:[%s,%s]",
					pair.ligand.sequence.toString(),
					formatBig(pair.complex.pfuncLowerBound),
					formatBig(pair.complex.pfuncUpperBoundSampled),
					formatBig(pair.ligand.pfuncLowerBound),
					formatBig(pair.ligand.pfuncUpperBoundSampled)
				);
				if (pair.ligand.sequence.isWildType()) {
					box.extraStyle = wildtypeBoxStyle;
				}

				log("%s\n%s", pair.complex, pair.ligand);
			}
			boxes.xaxis = boxes.makeXAxis();
			boxes.xaxis.tickFormat = "%.0f";
			boxes.yaxis = boxes.makeYAxis();
			boxes.yaxis.tickFormat = "%.0f";
			boxes.draw(svg);
			boxes.setBounds(svg, 10, 16);

			svg.finish().write(new File("pfuncs.svg"));
		}

		// sort sequences by decreasing binding score upper bounds
		pairs.sort(Comparator.comparing((SequenceInfoPair pair) ->
			MathTools.bigDivide(pair.complex.pfuncUpperBound, pair.ligand.pfuncLowerBound, decimalPrecision)
		).reversed());

		// plot bounds on binding score
		{
			SVG svg = new SVG();

			SVG.StyleClass wildtypeIntervalStyle = svg.makeStyleClass("wildtype-interval-style");
			wildtypeIntervalStyle.priority = 10;
			wildtypeIntervalStyle.setStrokeColor(0x66cc55);
			wildtypeIntervalStyle.setStrokeWidth(0.5);

			SVGPlot.Intervals intervals = new SVGPlot.Intervals();
			intervals.intervalWidth = 2.0;
			for (SequenceInfoPair pair : pairs) {

				// compute the binding score
				BigDecimal lowerBinding = MathTools.bigDivide(pair.complex.pfuncLowerBound, pair.ligand.pfuncUpperBound, decimalPrecision);
				BigDecimal upperBinding = MathTools.bigDivide(pair.complex.pfuncUpperBound, pair.ligand.pfuncLowerBound, decimalPrecision);

				SVGPlot.Intervals.Interval interval = intervals.addInterval(
					MathTools.log10p1(lowerBinding),
					MathTools.log10p1(upperBinding)
				);
				interval.id = String.format("%s: [%s,%s]",
					pair.ligand.sequence.toString(),
					formatBig(lowerBinding),
					formatBig(upperBinding)
				);
				if (pair.ligand.sequence.isWildType()) {
					interval.extraStyle = wildtypeIntervalStyle;
				}
			}
			intervals.axis = intervals.makeAxis();
			intervals.axis.tickFormat = "%.0f";
			intervals.draw(svg);
			intervals.setBounds(svg, 10, 16);

			svg.finish().write(new File("binding.svg"));
		}

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

		private double getPfuncUpperBoundFractionSampled() {
			return MathTools.bigDivide(pfuncUpperBoundSampled, pfuncUpperBound, decimalPrecision).doubleValue();
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
				log(buf, "\t\tpercent sampled:                 %.4f", getPfuncUpperBoundFractionSampled()*100.0);
			}
			if (pfuncLowerBound != null) {
				log(buf, "\tpfunc lower bound:             %s", formatBig(pfuncLowerBound));
				double logUncertainty = MathTools.log10p1(pfuncUpperBound) - MathTools.log10p1(pfuncLowerBound);
				log(buf, "\tpfunc log uncertainty:         %.2f", logUncertainty);
			}
			return buf.toString();
		}
	}

	private static ConfAStarTree makeAStar(EnergyMatrix emat, RCs rcs) {
		ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
			.setTraditional()
			/*.setMPLP(new ConfAStarTree.MPLPBuilder()
				.setUpdater(new EdgeUpdater())
				.setNumIterations(40)
			)
			*/
			.build();
		astar.setParallelism(Parallelism.makeCpu(6));
		return astar;
	}
}
