package edu.duke.cs.osprey.pruning;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.math.BigInteger;
import java.util.function.Consumer;

/**
 * Simple, fast implementation of DEE.
 * Doesn't support triples yet
 */
public class SimpleDEE {

	// TODO: implement parallelism?

	private static class Reporter {

		int numSingles = 0;
		int numPairs = 0;
		int numSinglesPruned = 0;
		int numPairsPruned = 0;
		int numSinglesPrunedThisRound = 0;
		int numPairSPrunedThisRound = 0;

		public Reporter(SimpleConfSpace confSpace) {
			numSingles = confSpace.getNumResConfs();
			numPairs = confSpace.getNumResConfPairs();
		}

		public void updateCounts(PruningMatrix pmat) {

			// count singles
			int newNumSinglesPruned = pmat.countPrunedRCs();
			numSinglesPrunedThisRound = newNumSinglesPruned - numSinglesPruned;
			numSinglesPruned = newNumSinglesPruned;

			// count pairs
			int newNumPairsPruned = pmat.countPrunedPairs();
			numPairSPrunedThisRound = newNumPairsPruned - numPairsPruned;
			numPairsPruned = newNumPairsPruned;
		}

		public void report(String prefix) {

			System.out.print(String.format("%30s: ", prefix));

			// report singles
			System.out.print(String.format("pruned %d singles,  %d/%d (%.1f%%) total",
				numSinglesPrunedThisRound,
				numSinglesPruned, numSingles, 100.0f*numSinglesPruned/numSingles
			));

			System.out.print(String.format("\n%30s: ", ""));

			// report pairs
			System.out.print(String.format("pruned %d pairs,  %d/%d (%.1f%%) total",
				numPairSPrunedThisRound,
				numPairsPruned, numPairs, 100.0f*numPairsPruned/numPairs
			));

			System.out.println();
		}

		public boolean prunedAnythingThisRound() {
			return numSinglesPrunedThisRound > 0 || numPairSPrunedThisRound > 0;
		}
	}

	/**
	 * Runs Dead-End Elimination (DEE) to remove tuples of RCs that will not appear in low-energy conformations.
	 */
	public static class Runner {

		private Double singlesThreshold = 100.0;
		private Double pairsThreshold = 100.0;
		private Double singlesGoldsteinDiffThreshold = null;
		private Double pairsGoldsteinDiffThreshold = null;
		private boolean typeDependent = false;
		private int numIterations = Integer.MAX_VALUE;
		private boolean showProgress = false;

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

		public Runner setGoldsteinDiffThreshold(Double val) {
			setSinglesGoldsteinDiffThreshold(val);
			setPairsGoldsteinDiffThreshold(val);
			return this;
		}

		public Runner setTypeDependent(boolean val) {
			typeDependent = val;
			return this;
		}

		public Runner setNumIterations(int val) {
			numIterations = val;
			return this;
		}

		public Runner setShowProgress(boolean val) {
			showProgress = val;
			return this;
		}

		public PruningMatrix run(SimpleConfSpace confSpace, EnergyMatrix emat) {

			if (showProgress) {
				System.out.println("Running DEE...");
			}

			Reporter reporter = new Reporter(confSpace);
			SimpleDEE dee = new SimpleDEE(confSpace, emat);

			Consumer<String> maybeReport = (prefix) -> {
				if (showProgress) {
					reporter.updateCounts(dee.pmat);
					reporter.report(prefix);
				}
			};

			// 1. threshold pruning
			if (singlesThreshold != null || pairsThreshold != null) {
				if (singlesThreshold != null) {
					dee.pruneSinglesByThreshold(singlesThreshold);
					maybeReport.accept("Threshold Singles");
				}
				if (pairsThreshold != null) {
					dee.prunePairsByThreshold(pairsThreshold);
					maybeReport.accept("Threshold Pairs");
				}
			}

			// 2. iterative DEE pruning
			if (singlesGoldsteinDiffThreshold != null || pairsGoldsteinDiffThreshold != null) {
				for (int i = 0; i<numIterations; i++) {

					if (showProgress) {
						System.out.println("DEE iteration " + (i+1) + "...");
					}

					int numPruned = reporter.numSinglesPruned + reporter.numPairsPruned;

					// 2.1 Goldstein criterion
					if (singlesGoldsteinDiffThreshold != null) {
						dee.pruneSinglesGoldstein(singlesGoldsteinDiffThreshold, typeDependent);
						maybeReport.accept("Goldstein Singles");
					}
					if (pairsGoldsteinDiffThreshold != null) {
						dee.prunePairsGoldstein(pairsGoldsteinDiffThreshold, typeDependent);
						maybeReport.accept("Goldstein Pairs");
					}

					// TODO: other pruning criteria?

					// stop if we didn't prune anything
					int numPrunedThisRound = reporter.numSinglesPruned + reporter.numPairsPruned - numPruned;
					if (numPrunedThisRound <= 0) {
						break;
					}
				}
			}

			// show total pruning results
			if (showProgress) {
				BigInteger allConfs = new RCs(confSpace).getNumConformations();
				BigInteger remainingConfs = new RCs(dee.pmat).getNumConformations();
				BigInteger prunedConfs = allConfs.subtract(remainingConfs);
				System.out.println("Conformations defined by conformation space:    " + String.format("%e", allConfs.floatValue()));
				System.out.println("Conformations remaining after pruning singles:  " + String.format("%e", remainingConfs.floatValue()));
				System.out.println("Percent conformations pruned by singles:        " + 100.0*prunedConfs.doubleValue()/allConfs.doubleValue());
			}

			return dee.pmat;
		}
	}

	public final SimpleConfSpace confSpace;
	public final EnergyMatrix emat;
	public final PruningMatrix pmat;

	public SimpleDEE(SimpleConfSpace confSpace, EnergyMatrix emat) {

		this.confSpace = confSpace;
		this.emat = emat;
		this.pmat = new PruningMatrix(confSpace);
	}

	public static enum SearcherResult {
		FoundIt,
		KeepSearching
	}

	public static interface SingleSearcher {
		public SearcherResult check(SimpleConfSpace.Position pos, SimpleConfSpace.ResidueConf rc);
	}

	public boolean findUnprunedSingle(SingleSearcher searcher) {
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			if (findUnprunedSingleAt(pos, searcher)) {
				return true;
			}
		}
		return false;
	}

	public boolean findUnprunedSingleAt(SimpleConfSpace.Position pos, SingleSearcher searcher) {
		for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {

			// skip pruned stuff
			if (pmat.isSinglePruned(pos.index, rc.index)) {
				continue;
			}

			if (searcher.check(pos, rc) == SearcherResult.FoundIt) {
				return true;
			}
		}
		return false;
	}

	public static interface SingleProcessor {
		public void process(SimpleConfSpace.Position pos, SimpleConfSpace.ResidueConf rc);

		public default SingleSearcher toSearcher() {
			return (pos, rc) -> {
				process(pos, rc);
				return SearcherResult.KeepSearching;
			};
		}
	}

	public void forEachUnprunedSingle(SingleProcessor processor) {
		findUnprunedSingle(processor.toSearcher());
	}

	public void forEachUnprunedSingleAt(SimpleConfSpace.Position pos, SingleProcessor processor) {
		findUnprunedSingleAt(pos, processor.toSearcher());
	}

	public static interface PairSearcher {
		public SearcherResult check(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf r2);
	}

	public boolean findUnprunedPair(PairSearcher searcher) {
		for (SimpleConfSpace.Position pos1 : confSpace.positions) {
			for (SimpleConfSpace.Position pos2 : confSpace.positions) {

				// only pairs, not cartesian product
				if (pos2.index >= pos1.index) {
					break;
				}

				if (findUnprunedPairAt(pos1, pos2, searcher)) {
					return true;
				}
			}
		}
		return false;
	}

	public boolean findUnprunedPairAt(SimpleConfSpace.Position pos1, SimpleConfSpace.Position pos2, PairSearcher searcher) {
		for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {

			// skip pruned stuff
			if (pmat.isSinglePruned(pos1.index, rc1.index)) {
				continue;
			}

			for (SimpleConfSpace.ResidueConf rc2 : pos2.resConfs) {

				// skip pruned stuff
				if (pmat.isPairPruned(pos1.index, rc1.index, pos2.index, rc2.index)) {
					continue;
				}

				if (searcher.check(pos1, rc1, pos2, rc2) == SearcherResult.FoundIt) {
					return true;
				}
			}
		}
		return false;
	}

	public static interface PairProcessor {
		public void process(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf r2);

		public default PairSearcher toSearcher() {
			return (pos1, rc1, pos2, rc2) -> {
				process(pos1, rc1, pos2, rc2);
				return SearcherResult.KeepSearching;
			};
		}
	}

	public void forEachUnprunedPair(PairProcessor processor) {
		findUnprunedPair(processor.toSearcher());
	}

	public void forEachUnprunedPairAt(SimpleConfSpace.Position pos1, SimpleConfSpace.Position pos2, PairProcessor processor) {
		findUnprunedPairAt(pos1, pos2, processor.toSearcher());
	}

	public static interface SinglePruner {
		public boolean shouldPrune(SimpleConfSpace.Position pos, SimpleConfSpace.ResidueConf rc);
	}

	public void pruneSingles(SinglePruner pruner) {
		forEachUnprunedSingle((pos, rc) -> {
			if (pruner.shouldPrune(pos, rc)) {
				pmat.pruneSingle(pos.index, rc.index);
			}
		});
	}

	public static interface PairPruner {
		public boolean shouldPrune(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf r2);
	}

	public void prunePairs(PairPruner pruner) {
		forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			if (pruner.shouldPrune(pos1, rc1, pos2, rc2)) {
				pmat.prunePair(pos1.index, rc1.index, pos2.index, rc2.index);
			}
		});
	}

	public void pruneSinglesByThreshold(double energyThreshold) {
		pruneSingles((pos, rc) -> {
			return emat.getOneBody(pos.index, rc.index) > energyThreshold;
		});
	}

	public void prunePairsByThreshold(double energyThreshold) {
		prunePairs((pos1, rc1, pos2, rc2) -> {
			return emat.getPairwise(pos1.index, rc1.index, pos2.index, rc2.index) > energyThreshold;
		});
	}

	public void pruneSinglesGoldstein(double energyDiffThreshold, boolean typeDependent) {
		pruneSingles((candidatePos, candidateRc) -> {

			// prune the candidate rc if we can find a competitor rc that has much lower energy
			return findUnprunedSingleAt(candidatePos, (competitorPos, competitorRc) -> {

				// don't compete against self
				if (competitorRc == candidateRc) {
					return SearcherResult.KeepSearching;
				}

				// skip unmatched types if needed
				if (typeDependent) {
					if (!candidateRc.template.name.equals(competitorRc.template.name)) {
						return SearcherResult.KeepSearching;
					}
				}

				// start with singles energy diff
				double energyDiffSum = 0
					+ emat.getOneBody(candidatePos.index, candidateRc.index)
					- emat.getOneBody(competitorPos.index, competitorRc.index);

				// sum over witness positions
				for (SimpleConfSpace.Position witnessPos : confSpace.positions) {

					// witness pos can't be candidate pos
					if (witnessPos == candidatePos) {
						continue;
					}

					// min over witness rcs
					double minEnergyDiff = Double.POSITIVE_INFINITY;
					for (SimpleConfSpace.ResidueConf witnessRc : witnessPos.resConfs) {

						// skip pruned witnesses
						if (pmat.isPairPruned(candidatePos.index, candidateRc.index, witnessPos.index, witnessRc.index)) {
							continue;
						}

						// compute the energy diff between the candidate and competitor, from the point of view of the witness
						double energyDiff = 0
							+ emat.getPairwise(candidatePos.index, candidateRc.index, witnessPos.index, witnessRc.index)
							- emat.getPairwise(competitorPos.index, competitorRc.index, witnessPos.index, witnessRc.index);
						minEnergyDiff = Math.min(minEnergyDiff, energyDiff);
					}

					energyDiffSum += minEnergyDiff;
					if (energyDiffSum == Double.POSITIVE_INFINITY) {
						break;
					}
				}

				if (energyDiffSum > energyDiffThreshold) {
					return SearcherResult.FoundIt;
				} else {
					return SearcherResult.KeepSearching;
				}
			});
		});
	}

	public void prunePairsGoldstein(double energyDiffThreshold, boolean typeDependent) {
		prunePairs((candidatePos1, candidateRc1, candidatePos2, candidateRc2) -> {

			// prune the candidate rc if we can find a competitor rc that has much lower energy
			return findUnprunedPairAt(candidatePos1, candidatePos2, (competitorPos1, competitorRc1, competitorPos2, competitorRc2) -> {

				// don't compete against self
				if (competitorRc1 == candidateRc1 && competitorRc2 == candidateRc2) {
					return SearcherResult.KeepSearching;
				}

				// skip unmatched types if needed
				if (typeDependent) {
					if (!candidateRc1.template.name.equals(competitorRc1.template.name)
						|| !candidateRc2.template.name.equals(competitorRc2.template.name)) {
						return SearcherResult.KeepSearching;
					}
				}

				// start with fragment energy diff
				double energyDiffSum = 0
					+ emat.getOneBody(candidatePos1.index, candidateRc1.index)
					+ emat.getOneBody(candidatePos2.index, candidateRc2.index)
					+ emat.getPairwise(candidatePos1.index, candidateRc1.index, candidatePos2.index, candidateRc2.index)
					- emat.getOneBody(competitorPos1.index, competitorRc1.index)
					- emat.getOneBody(competitorPos2.index, competitorRc2.index)
					- emat.getPairwise(competitorPos1.index, competitorRc1.index, competitorPos2.index, competitorRc2.index);

				// sum over witness positions
				for (SimpleConfSpace.Position witnessPos : confSpace.positions) {

					// witness pos can't be candidate pos
					if (witnessPos == candidatePos1 || witnessPos == candidatePos2) {
						continue;
					}

					// min over witness rcs
					double minEnergyDiff = Double.POSITIVE_INFINITY;
					for (SimpleConfSpace.ResidueConf witnessRc : witnessPos.resConfs) {

						// skip pruned witnesses
						if (pmat.isPairPruned(candidatePos1.index, candidateRc1.index, witnessPos.index, witnessRc.index)
							|| pmat.isPairPruned(candidatePos2.index, candidateRc2.index, witnessPos.index, witnessRc.index)) {
							continue;
						}

						// compute the energy diff between the candidate and competitor, from the point of view of the witness
						double energyDiff = 0
							+ emat.getPairwise(candidatePos1.index, candidateRc1.index, witnessPos.index, witnessRc.index)
							+ emat.getPairwise(candidatePos2.index, candidateRc2.index, witnessPos.index, witnessRc.index)
							- emat.getPairwise(competitorPos1.index, competitorRc1.index, witnessPos.index, witnessRc.index)
							- emat.getPairwise(competitorPos2.index, competitorRc2.index, witnessPos.index, witnessRc.index);
						minEnergyDiff = Math.min(minEnergyDiff, energyDiff);
					}
					energyDiffSum += minEnergyDiff;
				}

				if (energyDiffSum > energyDiffThreshold) {
					return SearcherResult.FoundIt;
				} else {
					return SearcherResult.KeepSearching;
				}
			});
		});
	}
}
