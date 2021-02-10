/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.pruning;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Progress;

import java.io.File;
import java.math.BigInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Consumer;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Simple, fast implementation of DEE.
 */
public class SimpleDEE {

	private static class Reporter {

		int numSingles;
		int numPairs;
		int numTriples;
		int numSinglesPruned = 0;
		int numPairsPruned = 0;
		int numTriplesPruned = 0;
		int numSinglesPrunedThisRound = 0;
		int numPairsPrunedThisRound = 0;
		int numTriplesPrunedThisRound = 0;

		public Reporter(SimpleConfSpace confSpace) {
			numSingles = confSpace.getNumResConfs();
			numPairs = confSpace.getNumResConfPairs();
			numTriples = confSpace.getNumResConfTriples();
		}

		public void updateCounts(PruningMatrix pmat) {

			// count singles
			int newNumSinglesPruned = pmat.countPrunedRCs();
			numSinglesPrunedThisRound = newNumSinglesPruned - numSinglesPruned;
			numSinglesPruned = newNumSinglesPruned;

			// count pairs
			int newNumPairsPruned = pmat.countPrunedPairs();
			numPairsPrunedThisRound = newNumPairsPruned - numPairsPruned;
			numPairsPruned = newNumPairsPruned;

			// count triples
			int newNumTriplesPruned = pmat.countPrunedTriples();
			numTriplesPrunedThisRound = newNumTriplesPruned - numTriplesPruned;
			numTriplesPruned = newNumTriplesPruned;
		}

		public void report(String prefix) {

			// report singles
			System.out.println(String.format("%30s: pruned %d singles,  %d/%d (%.1f%%) total",
				prefix,
				numSinglesPrunedThisRound,
				numSinglesPruned, numSingles, 100.0f*numSinglesPruned/numSingles
			));

			// report pairs
			System.out.println(String.format("%30s: pruned %d pairs,    %d/%d (%.1f%%) total",
				"",
				numPairsPrunedThisRound,
				numPairsPruned, numPairs, 100.0f*numPairsPruned/numPairs
			));

			// report triples
			System.out.println(String.format("%30s: pruned %d triples,  %d/%d (%.1f%%) total",
				"",
				numTriplesPrunedThisRound,
				numTriplesPruned, numTriples, 100.0f*numTriplesPruned/numTriples
			));
		}

		public boolean prunedAnythingThisRound() {
			return numSinglesPrunedThisRound > 0 || numPairsPrunedThisRound > 0 || numTriplesPrunedThisRound > 0;
		}
	}


	/**
	 * Reads a saved pruning matrix from disk, or throws an exception
	 */
	public static PruningMatrix read(SimpleConfSpace confSpace, File cacheFile) {
		return ObjectIO.readOrThrow(
			cacheFile,
			PruningMatrix.class,
			"pruning matrix",
			(pmat) -> pmat.matches(confSpace)
		);
	}

	/**
	 * Runs various pruning algorithms including
	 *    steric/threshold pruning
	 * 	  Dead-End Elimination (DEE)
	 *    Pruning of Local Unrealistic Geometries (PLUG)
	 *    Transitive pruning (ie, prune tuples implied by other pruned tuples)
	 * to remove tuples of RCs that will not appear in low-energy conformations.
	 */
	public static class Runner {

		private Double singlesThreshold = 100.0;
		private Double pairsThreshold = 100.0;
		private Double singlesGoldsteinDiffThreshold = null;
		private Double pairsGoldsteinDiffThreshold = null;
		private Double triplesGoldsteinDiffThreshold = null;
		private boolean typeDependent = false;
		private int numIterations = Integer.MAX_VALUE;
		private Double singlesPlugThreshold = null;
		private Double pairsPlugThreshold = null;
		private Double triplesPlugThreshold = null;
		private boolean singlesTransitivePruning = false;
		private boolean pairsTransitivePruning = false;
		private boolean triplesTransitivePruning = false;
		private boolean showProgress = false;
		private File cacheFile = null;
		private Parallelism parallelism = Parallelism.makeCpu(1);

		public Runner setSinglesThreshold(Double val) {
			singlesThreshold = val;
			return this;
		}

		public Runner setPairsThreshold(Double val) {
			pairsThreshold = val;
			return this;
		}

		public Runner setThreshold(Double val) {
			setSinglesThreshold(val);
			setPairsThreshold(val);
			return this;
		}

		public Runner setSinglesGoldsteinDiffThreshold(Double val) {
			singlesGoldsteinDiffThreshold = val;
			return this;
		}

		public Runner setPairsGoldsteinDiffThreshold(Double val) {
			pairsGoldsteinDiffThreshold = val;
			return this;
		}

		public Runner setTriplesGoldsteinDiffThreshold(Double val) {
			triplesGoldsteinDiffThreshold = val;
			return this;
		}

		public Runner setGoldsteinDiffThreshold(Double val) {
			setSinglesGoldsteinDiffThreshold(val);
			setPairsGoldsteinDiffThreshold(val);
			setTriplesGoldsteinDiffThreshold(val);
			return this;
		}

		public Runner setTypeDependent(boolean val) {
			typeDependent = val;

			// turn off threshold pruning since it's incompatible with type dependence
			singlesThreshold = null;
			pairsThreshold = null;

			return this;
		}

		public Runner setNumIterations(int val) {
			numIterations = val;
			return this;
		}

		public Runner setSinglesPlugThreshold(double val) {
			singlesPlugThreshold = val;
			return this;
		}

		public Runner setPairsPlugThreshold(double val) {
			pairsPlugThreshold = val;
			return this;
		}

		public Runner setTriplesPlugThreshold(double val) {
			triplesPlugThreshold = val;
			return this;
		}

		public Runner setPlugThreshold(double val) {
			setSinglesPlugThreshold(val);
			setPairsPlugThreshold(val);
			setTriplesPlugThreshold(val);
			return this;
		}

		public Runner setSinglesTransitivePruning(boolean val) {
			singlesTransitivePruning = val;
			return this;
		}

		public Runner setPairsTransitivePruning(boolean val) {
			pairsTransitivePruning = val;
			return this;
		}

		public Runner setTriplesTransitivePruning(boolean val) {
			triplesTransitivePruning = val;
			return this;
		}

		public Runner setTransitivePruning(boolean val) {
			setSinglesTransitivePruning(val);
			setPairsTransitivePruning(val);
			setTriplesTransitivePruning(val);
			return this;
		}

		public Runner setShowProgress(boolean val) {
			showProgress = val;
			return this;
		}

		public Runner setCacheFile(File val) {
			cacheFile = val;
			return this;
		}

		public Runner setParallelism(Parallelism val) {
			parallelism = val;
			return this;
		}

		public PruningMatrix run(SimpleConfSpace confSpace, EnergyMatrix emat) {

			// check the cache file first
			if (cacheFile != null) {
				return ObjectIO.readOrMake(
					cacheFile,
					PruningMatrix.class,
					"pruning matrix",
					(pmat) -> pmat.matches(confSpace),
					(context) -> reallyRun(confSpace, emat)
				);
			} else {
				return reallyRun(confSpace, emat);
			}
		}

		private PruningMatrix reallyRun(SimpleConfSpace confSpace, EnergyMatrix emat) {

			if (showProgress) {
				System.out.println("Running DEE...");
			}

			Reporter reporter = new Reporter(confSpace);
			PruningMatrix pmat = new PruningMatrix(confSpace);

			Consumer<String> maybeReport = (prefix) -> {
				if (showProgress) {
					reporter.updateCounts(pmat);
					reporter.report(prefix);
				}
			};

			// 1. threshold pruning
			if (singlesThreshold != null || pairsThreshold != null) {

				// can't use with type-dependent pruning
				if (typeDependent) {
					throw new IllegalArgumentException("threshold pruning can't be used with type dependence");
				}

				SimpleDEE dee = new SimpleDEE(confSpace, emat, pmat);
				if (singlesThreshold != null) {
					dee.pruneSinglesByThreshold(singlesThreshold);
					maybeReport.accept("Threshold Singles");
				}
				if (pairsThreshold != null) {
					dee.prunePairsByThreshold(pairsThreshold);
					maybeReport.accept("Threshold Pairs");
				}

				transitivePruning(confSpace, pmat, maybeReport);
			}

			// 2. iterative DEE pruning
			if (singlesGoldsteinDiffThreshold != null || pairsGoldsteinDiffThreshold != null || triplesGoldsteinDiffThreshold != null) {

				// choose the competitor RCs (prune with an interval of 0)
				if (showProgress) {
					System.out.println("Choosing competitor residue conformations...");
				}
				PruningMatrix competitors = new PruningMatrix(pmat);
				{
					SimpleDEE dee = new SimpleDEE(confSpace, emat, competitors);
					if (singlesGoldsteinDiffThreshold != null) {
						dee.pruneSinglesGoldstein(0, typeDependent);
					}
					if (pairsGoldsteinDiffThreshold != null) {
						dee.prunePairsGoldstein(0, typeDependent, parallelism);
					}
					if (triplesGoldsteinDiffThreshold != null) {
						dee.pruneTriplesGoldstein(0, typeDependent, parallelism);
						maybeReport.accept("Goldstein Triples");
					}
				}

				SimpleDEE dee = new SimpleDEE(confSpace, emat, pmat, competitors);

				for (int i=0; i<numIterations; i++) {

					if (showProgress) {
						System.out.println("DEE iteration " + (i+1) + "...");
					}

					int numPruned = reporter.numSinglesPruned + reporter.numPairsPruned + reporter.numTriplesPruned;

					// 3.1 Goldstein criterion
					if (singlesGoldsteinDiffThreshold != null) {
						dee.pruneSinglesGoldstein(singlesGoldsteinDiffThreshold, typeDependent);
						maybeReport.accept("Goldstein Singles");
					}
					if (pairsGoldsteinDiffThreshold != null) {
						dee.prunePairsGoldstein(pairsGoldsteinDiffThreshold, typeDependent, parallelism);
						maybeReport.accept("Goldstein Pairs");
					}
					if (triplesGoldsteinDiffThreshold != null) {
						dee.pruneTriplesGoldstein(triplesGoldsteinDiffThreshold, typeDependent, parallelism);
						maybeReport.accept("Goldstein Triples");
					}

					// TODO: other pruning criteria?

					// stop if we didn't prune anything
					int numPrunedThisRound = reporter.numSinglesPruned + reporter.numPairsPruned + reporter.numTriplesPruned - numPruned;
					if (numPrunedThisRound <= 0) {
						break;
					}
				}

				transitivePruning(confSpace, pmat, maybeReport);
			}

			// 3. pruning with PLUG
			if (singlesPlugThreshold != null || pairsPlugThreshold != null || triplesPlugThreshold != null) {
				if (showProgress) {
					System.out.println("Pruning with PLUG...");
				}
				try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {
					PLUG plug = new PLUG(confSpace);
					if (singlesPlugThreshold != null) {
						plug.pruneSingles(pmat, singlesPlugThreshold, tasks);
						maybeReport.accept("PLUG Singles");
					}
					if (pairsPlugThreshold != null) {
						plug.prunePairs(pmat, pairsPlugThreshold, tasks);
						maybeReport.accept("PLUG Pairs");
					}
					if (triplesPlugThreshold != null) {
						plug.pruneTriples(pmat, triplesPlugThreshold, tasks);
						maybeReport.accept("PLUG Triples");
					}
				}

				transitivePruning(confSpace, pmat, maybeReport);
			}

			// show total pruning results
			if (showProgress) {
				BigInteger allConfs = new RCs(confSpace).getNumConformations();
				BigInteger unprunedLower = pmat.calcUnprunedConfsLowerBound();
				BigInteger unprunedUpper = pmat.calcUnprunedConfsUpperBound();
				BigInteger prunedLower = allConfs.subtract(unprunedUpper);
				BigInteger prunedUpper = allConfs.subtract(unprunedLower);
				double percentPrunedLower = 100.0*prunedLower.doubleValue()/allConfs.doubleValue();
				double percentPrunedUpper = 100.0*prunedUpper.doubleValue()/allConfs.doubleValue();
				log("Conformations defined by conformation space:                           %s", formatBig(allConfs));
				log("Conformations pruned (by singles and pairs, bounds):                   [%s,%s]", formatBig(prunedLower), formatBig(prunedUpper));
				log("Conformations remaining after pruning (by singles and pairs, bounds):  [%s,%s]", formatBig(unprunedLower), formatBig(unprunedUpper));
				log("Percent conformations pruned (by singles and pairs, bounds):           [%.6f,%.6f]", percentPrunedLower, percentPrunedUpper);
			}

			return pmat;
		}

		private void transitivePruning(SimpleConfSpace confSpace, PruningMatrix pmat, Consumer<String> maybeReport) {
			if (singlesTransitivePruning || pairsTransitivePruning || triplesTransitivePruning) {
				if (showProgress) {
					System.out.println("Transitive Pruning...");
				}
				try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {
					TransitivePruner pruner = new TransitivePruner(confSpace);
					if (singlesTransitivePruning) {
						pruner.pruneSingles(pmat, tasks);
						maybeReport.accept("Transitive Singles");
					}
					if (pairsTransitivePruning) {
						pruner.prunePairs(pmat, tasks);
						maybeReport.accept("Transitive Pairs");
					}
					if (triplesTransitivePruning) {
						pruner.pruneTriples(pmat, tasks);
						maybeReport.accept("Transitive Triples");
					}
				}
			}
		}
	}

	public final SimpleConfSpace confSpace;
	public final EnergyMatrix emat;
	public final PruningMatrix pmat;
	public final PruningMatrix competitors;

	public SimpleDEE(SimpleConfSpace confSpace, EnergyMatrix emat, PruningMatrix pmat) {
		this(confSpace, emat, pmat, pmat);
	}

	public SimpleDEE(SimpleConfSpace confSpace, EnergyMatrix emat, PruningMatrix pmat, PruningMatrix competitors) {
		this.confSpace = confSpace;
		this.emat = emat;
		this.pmat = pmat;
		this.competitors = competitors;
	}

	private ResidueTemplate getTemplate(int pos, int rc) {
		return confSpace.positions.get(pos).resConfs.get(rc).template;
	}

	public void pruneSinglesByThreshold(double energyThreshold) {
		pmat.forEachUnprunedSingle((pos, rc) -> {
			if (emat.getOneBody(pos, rc) > energyThreshold) {
				pmat.pruneSingle(pos, rc);
			}
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	public void prunePairsByThreshold(double energyThreshold) {
		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			if (emat.getPairwise(pos1, rc1, pos2, rc2) > energyThreshold) {
				pmat.prunePair(pos1, rc1, pos2, rc2);
			}
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	public void pruneSinglesGoldstein(double energyDiffThreshold, boolean typeDependent) {

		// singles are so fast, we don't need to bother with parallelism and progress (right?)

		pmat.forEachUnprunedSingle((candidatePos, candidateRc) -> {

			// is there a competitor rc that has much lower energy?
			PruningMatrix.IteratorCommand result = competitors.forEachUnprunedSingleAt(candidatePos, (competitorPos, competitorRc) -> {

				// don't compete against self
				if (competitorRc == candidateRc) {
					return PruningMatrix.IteratorCommand.Continue;
				}

				// skip unmatched types if needed
				if (typeDependent) {
					if (!getTemplate(candidatePos, candidateRc).name.equals(getTemplate(competitorPos, competitorRc).name)) {
						return PruningMatrix.IteratorCommand.Continue;
					}
				}

				// start with singles energy diff
				double energyDiffSum = 0
					+ emat.getOneBody(candidatePos, candidateRc)
					- emat.getOneBody(competitorPos, competitorRc);

				// sum over witness positions
				for (int witnessPos=0; witnessPos<confSpace.positions.size(); witnessPos++) {

					// witness pos can't be candidate pos
					if (witnessPos == candidatePos) {
						continue;
					}

					// min over witness rcs
					double minEnergyDiff = Double.POSITIVE_INFINITY;
					int numWitnessRCs = confSpace.positions.get(witnessPos).resConfs.size();
					for (int witnessRc=0; witnessRc<numWitnessRCs; witnessRc++) {

						// skip pruned witnesses
						if (pmat.isPairPruned(candidatePos, candidateRc, witnessPos, witnessRc)) {
							continue;
						}

						// compute the energy diff between the candidate and competitor, from the point of view of the witness
						double energyDiff = 0
							+ emat.getPairwise(candidatePos, candidateRc, witnessPos, witnessRc)
							- emat.getPairwise(competitorPos, competitorRc, witnessPos, witnessRc);
						minEnergyDiff = Math.min(minEnergyDiff, energyDiff);
					}

					energyDiffSum += minEnergyDiff;
					if (energyDiffSum == Double.POSITIVE_INFINITY) {
						break;
					}
				}

				// if we found a suitable competitor, stop searching
				if (energyDiffSum > energyDiffThreshold) {
					return PruningMatrix.IteratorCommand.Break;
				} else {
					return PruningMatrix.IteratorCommand.Continue;
				}
			});

			// if the iteration terminated early (ie, we found a suitable competitor), then prune the candidate
			if (result == PruningMatrix.IteratorCommand.Break) {
				pmat.pruneSingle(candidatePos, candidateRc);
			}

			// always go onto to check the next candidate
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	public void prunePairsGoldstein(double energyDiffThreshold, boolean typeDependent) {
		prunePairsGoldstein(energyDiffThreshold, typeDependent, Parallelism.makeCpu(1));
	}

	public void prunePairsGoldstein(double energyDiffThreshold, boolean typeDependent, Parallelism parallelism) {

		// this one can take quite a while, so track progress
		AtomicLong numPairs = new AtomicLong(0);
		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			numPairs.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numPairs.get());

		try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {

			pmat.forEachUnprunedPair((candidatePos1, candidateRc1, candidatePos2, candidateRc2) -> {

				tasks.submit(
					() -> {
						// can we find a competitor rc that has much lower energy?
						return competitors.forEachUnprunedPairAt(candidatePos1, candidatePos2, (competitorPos1, competitorRc1, competitorPos2, competitorRc2) -> {

							// don't compete against self
							if (competitorRc1 == candidateRc1 && competitorRc2 == candidateRc2) {
								return PruningMatrix.IteratorCommand.Continue;
							}

							// skip unmatched types if needed
							if (typeDependent) {
								if (!getTemplate(candidatePos1, candidateRc1).name.equals(getTemplate(competitorPos1, competitorRc1).name)
									|| !getTemplate(candidatePos2, candidateRc2).name.equals(getTemplate(competitorPos2, competitorRc2).name)) {
									return PruningMatrix.IteratorCommand.Continue;
								}
							}

							// start with fragment energy diff
							double energyDiffSum = 0
								+ emat.getOneBody(candidatePos1, candidateRc1)
								+ emat.getOneBody(candidatePos2, candidateRc2)
								+ emat.getPairwise(candidatePos1, candidateRc1, candidatePos2, candidateRc2)
								- emat.getOneBody(competitorPos1, competitorRc1)
								- emat.getOneBody(competitorPos2, competitorRc2)
								- emat.getPairwise(competitorPos1, competitorRc1, competitorPos2, competitorRc2);

							// sum over witness positions
							for (int witnessPos=0; witnessPos<confSpace.positions.size(); witnessPos++) {

								// witness pos can't be candidate pos
								if (witnessPos == candidatePos1 || witnessPos == candidatePos2) {
									continue;
								}

								// min over witness rcs
								double minEnergyDiff = Double.POSITIVE_INFINITY;
								int numWitnessRCs = confSpace.positions.get(witnessPos).resConfs.size();
								for (int witnessRc=0; witnessRc<numWitnessRCs; witnessRc++) {

									// skip pruned witnesses
									if (pmat.isPairPruned(candidatePos1, candidateRc1, witnessPos, witnessRc)
										|| pmat.isPairPruned(candidatePos2, candidateRc2, witnessPos, witnessRc)) {
										continue;
									}

									// compute the energy diff between the candidate and competitor, from the point of view of the witness
									double energyDiff = 0
										+ emat.getPairwise(candidatePos1, candidateRc1, witnessPos, witnessRc)
										+ emat.getPairwise(candidatePos2, candidateRc2, witnessPos, witnessRc)
										- emat.getPairwise(competitorPos1, competitorRc1, witnessPos, witnessRc)
										- emat.getPairwise(competitorPos2, competitorRc2, witnessPos, witnessRc);
									minEnergyDiff = Math.min(minEnergyDiff, energyDiff);
								}
								energyDiffSum += minEnergyDiff;
								if (energyDiffSum == Double.POSITIVE_INFINITY) {
									break;
								}
							}

							if (energyDiffSum > energyDiffThreshold) {
								return PruningMatrix.IteratorCommand.Break;
							} else {
								return PruningMatrix.IteratorCommand.Continue;
							}
						});
					},
					(result) -> {

						// if the iteration terminated early (ie, we found a suitable competitor), then prune the candidate
						if (result == PruningMatrix.IteratorCommand.Break) {
							pmat.prunePair(candidatePos1, candidateRc1, candidatePos2, candidateRc2);
						}

						progress.incrementProgress();
					}
				);

				// always go onto to check the next candidate
				return PruningMatrix.IteratorCommand.Continue;
			});
		}
	}

	public void pruneTriplesGoldstein(double energyDiffThreshold, boolean typeDependent) {
		pruneTriplesGoldstein(energyDiffThreshold, typeDependent, Parallelism.makeCpu(1));
	}

	public void pruneTriplesGoldstein(double energyDiffThreshold, boolean typeDependent, Parallelism parallelism) {

		// this one can take quite a while, so track progress
		AtomicLong numTriples = new AtomicLong(0);
		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			numTriples.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numTriples.get());

		try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {

			pmat.forEachUnprunedTriple((candidatePos1, candidateRc1, candidatePos2, candidateRc2, candidatePos3, candidateRc3) -> {

				tasks.submit(
					() -> {
						// can we find a competitor rc that has much lower energy?
						return competitors.forEachUnprunedTripleAt(candidatePos1, candidatePos2, candidatePos3, (competitorPos1, competitorRc1, competitorPos2, competitorRc2, competitorPos3, competitorRc3) -> {

							// don't compete against self
							if (competitorRc1 == candidateRc1 && competitorRc2 == candidateRc2 && competitorRc3 == candidateRc3) {
								return PruningMatrix.IteratorCommand.Continue;
							}

							// skip unmatched types if needed
							if (typeDependent) {
								if (!getTemplate(candidatePos1, candidateRc1).name.equals(getTemplate(competitorPos1, competitorRc1).name)
									|| !getTemplate(candidatePos2, candidateRc2).name.equals(getTemplate(competitorPos2, competitorRc2).name)
									|| !getTemplate(candidatePos3, candidateRc3).name.equals(getTemplate(competitorPos3, competitorRc3).name)) {
									return PruningMatrix.IteratorCommand.Continue;
								}
							}

							// start with fragment energy diff
							double energyDiffSum = 0
								+ emat.getOneBody(candidatePos1, candidateRc1)
								+ emat.getOneBody(candidatePos2, candidateRc2)
								+ emat.getOneBody(candidatePos3, candidateRc3)
								+ emat.getPairwise(candidatePos1, candidateRc1, candidatePos2, candidateRc2)
								+ emat.getPairwise(candidatePos1, candidateRc1, candidatePos3, candidateRc3)
								+ emat.getPairwise(candidatePos2, candidateRc2, candidatePos3, candidateRc3)
								- emat.getOneBody(competitorPos1, competitorRc1)
								- emat.getOneBody(competitorPos2, competitorRc2)
								- emat.getOneBody(competitorPos3, competitorRc3)
								- emat.getPairwise(competitorPos1, competitorRc1, competitorPos2, competitorRc2)
								- emat.getPairwise(competitorPos1, competitorRc1, competitorPos3, competitorRc3)
								- emat.getPairwise(competitorPos2, competitorRc2, competitorPos3, competitorRc3);

							// sum over witness positions
							for (int witnessPos=0; witnessPos<confSpace.positions.size(); witnessPos++) {

								// witness pos can't be candidate pos
								if (witnessPos == candidatePos1 || witnessPos == candidatePos2 || witnessPos == candidatePos3) {
									continue;
								}

								// min over witness rcs
								double minEnergyDiff = Double.POSITIVE_INFINITY;
								int numWitnessRCs = confSpace.positions.get(witnessPos).resConfs.size();
								for (int witnessRc=0; witnessRc<numWitnessRCs; witnessRc++) {

									// skip pruned witnesses
									if (pmat.isPairPruned(candidatePos1, candidateRc1, witnessPos, witnessRc)
										|| pmat.isPairPruned(candidatePos2, candidateRc2, witnessPos, witnessRc)
										|| pmat.isPairPruned(candidatePos3, candidateRc3, witnessPos, witnessRc)) {
										continue;
									}

									// compute the energy diff between the candidate and competitor, from the point of view of the witness
									double energyDiff = 0
										+ emat.getPairwise(candidatePos1, candidateRc1, witnessPos, witnessRc)
										+ emat.getPairwise(candidatePos2, candidateRc2, witnessPos, witnessRc)
										+ emat.getPairwise(candidatePos3, candidateRc3, witnessPos, witnessRc)
										- emat.getPairwise(competitorPos1, competitorRc1, witnessPos, witnessRc)
										- emat.getPairwise(competitorPos2, competitorRc2, witnessPos, witnessRc)
										- emat.getPairwise(competitorPos3, competitorRc3, witnessPos, witnessRc);
									minEnergyDiff = Math.min(minEnergyDiff, energyDiff);
								}
								energyDiffSum += minEnergyDiff;
								if (energyDiffSum == Double.POSITIVE_INFINITY) {
									break;
								}
							}

							if (energyDiffSum > energyDiffThreshold) {
								return PruningMatrix.IteratorCommand.Break;
							} else {
								return PruningMatrix.IteratorCommand.Continue;
							}
						});
					},
					(result) -> {

						// if the iteration terminated early (ie, we found a suitable competitor), then prune the candidate
						if (result == PruningMatrix.IteratorCommand.Break) {
							pmat.pruneTriple(candidatePos1, candidateRc1, candidatePos2, candidateRc2, candidatePos3, candidateRc3);
						}

						progress.incrementProgress();
					}
				);

				// always go onto to check the next candidate
				return PruningMatrix.IteratorCommand.Continue;
			});
		}
	}
}
