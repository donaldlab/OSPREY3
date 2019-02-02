package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.*;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


public class SofeaLab {

	public static void main(String[] args) {

		ForcefieldParams ffparams = new ForcefieldParams();
		boolean recalc = false;

		// use the new templates, cuz why not
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.clearTemplateCoords()
			.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
			.build();

		// define design flexibility [68,73]
		Map<String,List<String>> designFlex = new HashMap<>();
		// unavoidable clash at A68. don't use ARG, or sub something smaller
		//designFlex.put("A68", Arrays.asList(Strand.WildType /* arg=34 */));
		//designFlex.put("A69", Arrays.asList(Strand.WildType /* ser=18 *//*, "THR", "LEU", "ILE", "VAL", "ALA", "GLY", "CYS"*/));
		designFlex.put("A70", Arrays.asList(Strand.WildType /* gly=1 */, "ALA", "VAL", "LEU", "ILE", "CYS"));
		//designFlex.put("A71", Arrays.asList(Strand.WildType /* lys=27 */));
		//designFlex.put("A72", Arrays.asList(Strand.WildType /* gln=9 */));
		designFlex.put("A73", Arrays.asList(Strand.WildType /* leu=5 */));

		// define target flexibility [5,10]
		List<String> targetFlex = Arrays.asList(
			//"A5", // lys=27
			//"A6", // hie=8
			//"A7", // tyr=8
			//"A8", // gln=9
			"A9", // phe=4
			"A10" // asn=7
		);

		// build strands
		Molecule pdb = PDBIO.readResource("/1CC8.ss.pdb");
		Strand design = new Strand.Builder(pdb)
			.setTemplateLibrary(templateLib)
			.setResidues("A68", "A73")
			.build();
		for (Map.Entry<String,List<String>> entry : designFlex.entrySet()) {
			design.flexibility.get(entry.getKey())
				.setLibraryRotamers(entry.getValue())
				.addWildTypeRotamers()
				.setContinuous();
		}
		Strand target = new Strand.Builder(pdb)
			.setTemplateLibrary(templateLib)
			.setResidues("A2", "A67")
			.build();
		for (String resNum : targetFlex) {
			target.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		// make a multi-state conf space
		Function<List<Strand>,SimpleConfSpace> makeConfSpace = (strands) ->
			new SimpleConfSpace.Builder().addStrands(strands)
				.setShellDistance(6.0)
				.build();
		MultiStateConfSpace confSpace = new MultiStateConfSpace
			.Builder("complex", makeConfSpace.apply(Arrays.asList(design, target)))
			.addMutableState("design", makeConfSpace.apply(Arrays.asList(design)))
			.addUnmutableState("target", makeConfSpace.apply(Arrays.asList(target)))
			.build();

		log("seq space: %s", confSpace.seqSpace);

		BiFunction<SimpleConfSpace,EnergyCalculator,ConfEnergyCalculator> makeConfEcalc = (simpleConfSpace, ecalc) ->
			new ConfEnergyCalculator.Builder(simpleConfSpace, ecalc)
				.setEnergyPartition(EnergyPartition.AllOnPairs) // use the tighter lower bounds
				.build();

		try (EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(2))
			.build()) {

			EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
				.setIsMinimizing(false)
				.build();

			// try a new minimizing kind of SOFEA
			Sofea sofea = new Sofea.Builder(confSpace)
				.configEachState(state -> {

					File ematLowerFile = new File(String.format("sofea.%s.lower.emat", state.name));
					File ematUpperFile = new File(String.format("sofea.%s.upper.emat", state.name));
					File confdbFile = new File(String.format("sofea.%s.confdb", state.name));
					if (recalc) {
						ematLowerFile.delete();
						ematUpperFile.delete();
						confdbFile.delete();
					}

					ConfEnergyCalculator minimizingConfEcalc = makeConfEcalc.apply(state.confSpace, minimizingEcalc);
					EnergyMatrix ematLower = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc)
						.setCacheFile(ematLowerFile)
						.build()
						.calcEnergyMatrix();

					ConfEnergyCalculator rigidConfEcalc = makeConfEcalc.apply(state.confSpace, rigidEcalc);
					EnergyMatrix ematUpper = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
						.setCacheFile(ematUpperFile)
						.build()
						.calcEnergyMatrix();

					return new Sofea.StateConfig(
						ematLower,
						ematUpper,
						minimizingConfEcalc,
						confdbFile
					);
				})
				.build();

			/* TEMP: brute force the free energies for all sequences
			for (Sequence seq : confSpace.seqSpace.getSequences()) {
				log("%s", seq);
				for (MultiStateConfSpace.State state : confSpace.states) {
					BigDecimal z = sofea.calcZSum(seq, state);
					double g = sofea.bcalc.freeEnergyPrecise(z);
					log("\t%10s   z=%s  g=%.4f", state.name, Log.formatBigLn(z), g);
				}
			}
			*/

			sofea.init(true);

			MultiStateConfSpace.LMFE lmfe = confSpace.lmfe()
				.addPositive("complex")
				.addNegative("design")
				.addNegative("target")
				.build();

			// TEMP: use the usual affinity optimization objective function
			//MinLMFE criterion = new MinLMFE(lmfe, 1);

			SequenceLMFE criterion = new SequenceLMFE(confSpace.seqSpace.makeWildTypeSequence(), lmfe, 0.5);

			// do the design!
			sofea.refine(criterion);
		}
	}

	static class SofeaTest {

		public final MathContext mathContext = new MathContext(64, RoundingMode.HALF_UP);
		public final BoltzmannCalculator bcalc = new BoltzmannCalculator(mathContext);

		class StateInfo {

			final MultiStateConfSpace.State state;
			final EnergyMatrix ematLower;
			final EnergyMatrix ematUpper;
			final ZMatrix zmatLower;
			final ZMatrix zmatUpper;

			public StateInfo(MultiStateConfSpace.State state, EnergyMatrix ematLower, EnergyMatrix ematUpper) {

				this.state = state;
				this.ematLower = ematLower;
				this.ematUpper = ematUpper;

				zmatLower = new ZMatrix(state.confSpace);
				zmatLower.set(ematUpper, bcalc);
				zmatUpper = new ZMatrix(state.confSpace);
				zmatUpper.set(ematLower, bcalc);
			}

			ConfIndex makeConfIndex() {
				ConfIndex index = new ConfIndex(state.confSpace.positions.size());
				index.updateUndefined();
				return index;
			}

			BigMath bigMath() {
				return new BigMath(mathContext);
			}

			BigDecimalBounds multBounds(BigDecimalBounds a, BigDecimalBounds b) {
				return new BigDecimalBounds(
					bigMath().set(a.lower).mult(b.lower).get(),
					bigMath().set(a.upper).mult(b.upper).get()
				);
			}

			BigDecimalBounds multBounds(BigDecimalBounds a, BigDecimal b) {
				return new BigDecimalBounds(
					bigMath().set(a.lower).mult(b).get(),
					bigMath().set(a.upper).mult(b).get()
				);
			}

			String bruteForceAstar(ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable, RCs rcs) {

				BigMath sum = bigMath().set(0.0);
				BigMath min = bigMath();
				BigMath max = bigMath();

				ConfAStarTree astar = new ConfAStarTree.Builder(ematLower, rcs)
					.setTraditional()
					.build();
				while (true) {

					ConfSearch.ScoredConf conf = astar.nextConf();
					if (conf == null) {
						break;
					}

					confEcalc.calcEnergyAsync(conf, confTable, (epmol) -> {
						BigDecimal z = bcalc.calcPrecise(epmol.getEnergy());
						sum.add(z);
						min.minOrSet(z);
						max.maxOrSet(z);
					});
				}

				confEcalc.ecalc.tasks.waitForFinish();

				return String.format("Zsum=%s  Zrange=[%s,%s]",
					Log.formatBigLn(sum.get()),
					Log.formatBigLn(min.get()),
					Log.formatBigLn(max.get())
				);
			}

			String bruteForce(ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable, RCs rcs) {
				return bruteForce(confEcalc, confTable, rcs, makeConfIndex());
			}

			String bruteForce(ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable, RCs rcs, ConfIndex index) {

				final BigMath sum = bigMath().set(0.0);
				final BigMath min = bigMath().set(MathTools.BigPositiveInfinity);
				final BigMath max = bigMath().set(MathTools.BigNegativeInfinity);

				final Runnable[] recurse = { null };
				recurse[0] = () -> {

					// are we fully assigned?
					if (index.isFullyDefined()) {

						// compute the conf energy
						RCTuple frag = new RCTuple(Conf.make(index));
						double e = confEcalc.calcEnergy(frag, confTable);
						BigDecimal z = bcalc.calcPrecise(e);
						sum.add(z);
						min.atMost(z);
						max.atLeast(z);

					} else {

						// otherwise, recurse
						int pos = index.numDefined;
						for (int rc : rcs.get(pos)) {
							index.assignInPlace(pos, rc);
							recurse[0].run();
							index.unassignInPlace(pos);
						}
					}
				};
				recurse[0].run();

				return String.format("Zsum=%s  Zrange=[%s,%s]",
					Log.formatBigLn(sum.get()),
					Log.formatBigLn(min.get()),
					Log.formatBigLn(max.get())
				);
			}

			void go(ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable, RCs rcs) {

				ConfIndex index = makeConfIndex();

				Runnable dump = () -> {

					BigDecimalBounds zPathHeadBounds = calcZPathHeadBounds(index);
					log("\tzPathHeadBounds = %s",
						Log.formatBigLn(zPathHeadBounds)
					);

					BigDecimalBounds zPathTailBounds = calcZPathTailBounds(index, rcs);
					BigDecimalBounds zPathBounds = multBounds(zPathHeadBounds, zPathTailBounds);
					log("\tzPathTailBounds = %s  zPathBounds = %s",
						Log.formatBigLn(zPathTailBounds),
						Log.formatBigLn(zPathBounds)
					);

					BigInteger leafNodes = countLeafNodes(index, rcs);
					BigDecimalBounds zSumTailBounds = new BigDecimalBounds(
						calcZSumTailLower(index, rcs),
						bigMath().set(zPathTailBounds.upper).mult(leafNodes).get()
					);
					BigDecimalBounds zSumBounds = multBounds(zPathHeadBounds, zSumTailBounds);
					log("\tzSumTailBounds = %s  zSumBounds = %s  leafNodes = %s",
						Log.formatBigLn(zSumTailBounds),
						Log.formatBigLn(zSumBounds),
						Log.formatBig(leafNodes)
					);

					if (index.isFullyDefined()) {

						BigDecimal zPath = calcZPath(index, confEcalc, confTable);
						log("\tzPath = %s",
							Log.formatBigLn(zPath)
						);

					} else {

						BigDecimalBounds zPathHeadBoundsTighter = calcZPathHeadBoundsTighter(index, confEcalc, confTable);
						log("\tzPathHeadBoundsTighter = %s",
							Log.formatBigLn(zPathHeadBoundsTighter)
						);

						BigDecimalBounds zPathBoundsTighter = multBounds(zPathHeadBoundsTighter, zPathTailBounds);
						log("\tzPathBoundsTighter = %s",
							Log.formatBigLn(zPathBoundsTighter)
						);

						BigDecimalBounds zSumBoundsTighter = multBounds(zPathHeadBoundsTighter, zSumTailBounds);
						log("\tzSumBoundsTighter = %s",
							Log.formatBigLn(zSumBoundsTighter)
						);

						zSumBoundsTighter.lower = calcZSumLowerTighter(index, rcs, confEcalc, confTable);
						log("\tzSumBoundsTighter = %s",
							Log.formatBigLn(zSumBoundsTighter)
						);
					}
				};

				// start at the root
				log("\nroot  %s", bruteForce(confEcalc, confTable, rcs, index));
				dump.run();

				int[] bestConf = {
					1, // max 1, best 1
					4, // max 5, best 5
					4, // max 4, best 4
					7 // max 7, best 7
				};

				for (int i=0; i<4; i++) {

					// make assignment at pos i
					index.assignInPlace(i, bestConf[i]);
					log("\nassignment %d  %s", i+1, bruteForce(confEcalc, confTable, rcs, index));

					dump.run();
				}
			}

			BigDecimalBounds calcZPathHeadBounds(ConfIndex index) {
				return new BigDecimalBounds(
					calcZPathHeadBound(zmatLower, index),
					calcZPathHeadBound(zmatUpper, index)
				);
			}

			BigDecimal calcZPathHeadBound(ZMatrix zmat, ConfIndex index) {

				BigMath z = bigMath().set(1.0);

				for (int i1=0; i1<index.numDefined; i1++) {
					int pos1 = index.definedPos[i1];
					int rc1 = index.definedRCs[i1];

					z.mult(zmat.getOneBody(pos1, rc1));

					for (int i2=0; i2<i1; i2++) {
						int pos2 = index.definedPos[i2];
						int rc2 = index.definedRCs[i2];

						z.mult(zmat.getPairwise(pos1, rc1, pos2, rc2));
					}
				}

				return z.get();
			}

			BigDecimal calcZPathNodeBound(ZMatrix zmat, ConfIndex index, int pos1, int rc1) {

				BigMath z = bigMath().set(zmat.getOneBody(pos1, rc1));

				for (int i=0; i<index.numDefined; i++) {
					int pos2 = index.definedPos[i];
					int rc2 = index.definedRCs[i];

					z.mult(zmat.getPairwise(pos1, rc1, pos2, rc2));
				}

				return z.get();
			}

			BigDecimalBounds calcZPathTailBounds(ConfIndex index, RCs rcs) {
				return new BigDecimalBounds(
					calcZPathTailBound(zmatLower, MathTools.Optimizer.Minimize, index, rcs),
					calcZPathTailBound(zmatUpper, MathTools.Optimizer.Maximize, index, rcs)
				);
			}

			BigDecimal calcZPathTailBound(ZMatrix zmat, MathTools.Optimizer opt, ConfIndex index, RCs rcs) {

				// this is the usual A* heuristic

				BigMath z = bigMath().set(1.0);

				// for each undefined position
				for (int i=0; i<index.numUndefined; i++) {
					int pos1 = index.undefinedPos[i];

					// optimize over possible assignments to pos1
					BigDecimal zpos1 = opt.initBigDecimal();
					for (int rc1 : rcs.get(pos1)) {

						BigMath zrc1 = bigMath()
							.set(zmat.getOneBody(pos1, rc1));

						// interactions with defined residues
						for (int j=0; j<index.numDefined; j++) {
							int pos2 = index.definedPos[j];
							int rc2 = index.definedRCs[j];

							zrc1.mult(zmat.getPairwise(pos1, rc1, pos2, rc2));
						}

						// interactions with undefined residues
						for (int j=0; j<i; j++) {
							int pos2 = index.undefinedPos[j];

							// optimize over possible assignments to pos2
							BigDecimal optzrc2 = opt.initBigDecimal();
							for (int rc2 : rcs.get(pos2)) {

								// pair with pos2
								BigDecimal zrc2 = zmat.getPairwise(pos1, rc1, pos2, rc2);

								optzrc2 = opt.opt(optzrc2, zrc2);
							}

							zrc1.mult(optzrc2);
						}

						zpos1 = opt.opt(zpos1, zrc1.get());
					}

					assert (MathTools.isFinite(zpos1));
					z.mult(zpos1);
				}

				return z.get();
			}

			BigInteger countLeafNodes(ConfIndex index, RCs rcs) {

				BigInteger count = BigInteger.ONE;

				for (int i=0; i<index.numUndefined; i++) {
					int pos = index.undefinedPos[i];
					count = count.multiply(BigInteger.valueOf(rcs.getNum(pos)));
				}

				return count;
			}

			void greedilyAssignConf(ConfIndex index, RCs rcs, ZMatrix zmat) {

				// get to any leaf node and compute the subtree part of its zPath

				int numUnassigned = index.numUndefined;

				// assign each unassigned position greedily
				for (int i=0; i<numUnassigned; i++) {

					int pos = index.numDefined;

					// find the RC with the biggest zPathComponent
					int bestRC = -1;
					BigDecimal bestZPathNode = MathTools.BigNegativeInfinity;
					for (int rc : rcs.get(pos)) {

						BigDecimal zPathNode = calcZPathNodeBound(zmat, index, pos, rc);
						if (MathTools.isGreaterThan(zPathNode, bestZPathNode)) {
							bestZPathNode = zPathNode;
							bestRC = rc;
						}
					}

					// make the assignment
					index.assignInPlace(pos, bestRC);
				}
			}

			void unassignConf(ConfIndex index, int numAssigned) {

				// undo all the assignments in the reverse order they were assigned
				for (int pos=index.numPos-1; pos>=numAssigned; pos--) {
					index.unassignInPlace(pos);
				}
			}

			BigDecimal calcZSumTailLower(ConfIndex index, RCs rcs) {

				// get to any leaf node and compute its zPathTailLower
				int numAssigned = index.numDefined;
				greedilyAssignConf(index, rcs, zmatLower);
				BigMath m = bigMath().set(1.0);
				for (int pos=index.numPos-1; pos>=numAssigned; pos--) {
					int rc = index.definedRCs[index.findDefined(pos)];
					index.unassignInPlace(pos);
					m.mult(calcZPathNodeBound(zmatLower, index, pos, rc));
				}
				return m.get();
			}

			BigDecimalBounds calcZPathHeadBoundsTighter(ConfIndex index, ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable) {
				return new BigDecimalBounds(
					calcZPathHeadLowerTighter(index, confEcalc, confTable),
					calcZPathHeadBound(zmatUpper, index)
				);
			}

			BigDecimal calcZPathHeadLowerTighter(ConfIndex index, ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable) {
				RCTuple tuple = new RCTuple(index);
				ResidueInteractions inters = confEcalc.makeTupleInters(tuple);
				double e = confEcalc.calcEnergy(tuple, inters, confTable);
				return bcalc.calcPrecise(e);
			}

			BigDecimal calcZPath(ConfIndex index, ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable) {

				if (!index.isFullyDefined()) {
					throw new IllegalArgumentException("not a full conf");
				}

				double e = confEcalc.calcEnergy(new RCTuple(index), confTable);
				return bcalc.calcPrecise(e);
			}

			BigDecimal calcZSumLowerTighter(ConfIndex index, RCs rcs, ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable) {

				// get to any leaf node and compute its zPath
				int numAssigned = index.numDefined;
				greedilyAssignConf(index, rcs, zmatUpper);
				BigDecimal zPath = calcZPath(index, confEcalc, confTable);
				unassignConf(index, numAssigned);
				return zPath;
			}
		}
	}
}
