package edu.duke.cs.osprey.pruning;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.restypes.ResidueTemplate;

import java.math.BigInteger;
import java.util.function.Consumer;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;

/**
 * Simple, fast implementation of DEE.
 */
public class SimpleDEE {

	// TODO: implement parallelism?

	private static class Reporter {

		int numSingles = 0;
		int numPairs = 0;
		int numTriples = 0;
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
	 * Runs Dead-End Elimination (DEE) to remove tuples of RCs that will not appear in low-energy conformations.
	 */
	public static class Runner {

		private Double singlesThreshold = 100.0;
		private Double pairsThreshold = 100.0;
		private Double singlesGoldsteinDiffThreshold = null;
		private Double pairsGoldsteinDiffThreshold = null;
		private Double triplesGoldsteinDiffThreshold = null;
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

					int numPruned = reporter.numSinglesPruned + reporter.numPairsPruned + reporter.numTriplesPruned;

					// 2.1 Goldstein criterion
					if (singlesGoldsteinDiffThreshold != null) {
						dee.pruneSinglesGoldstein(singlesGoldsteinDiffThreshold, typeDependent);
						maybeReport.accept("Goldstein Singles");
					}
					if (pairsGoldsteinDiffThreshold != null) {
						dee.prunePairsGoldstein(pairsGoldsteinDiffThreshold, typeDependent);
						maybeReport.accept("Goldstein Pairs");
					}
					if (triplesGoldsteinDiffThreshold != null) {
						dee.pruneTriplesGoldstein(triplesGoldsteinDiffThreshold, typeDependent);
						maybeReport.accept("Goldstein Triples");
					}

					// TODO: other pruning criteria?

					// stop if we didn't prune anything
					int numPrunedThisRound = reporter.numSinglesPruned + reporter.numPairsPruned + reporter.numTriplesPruned - numPruned;
					if (numPrunedThisRound <= 0) {
						break;
					}
				}
			}

			// show total pruning results
			if (showProgress) {
				BigInteger allConfs = new RCs(confSpace).getNumConformations();
				BigInteger unprunedLower = dee.pmat.calcUnprunedConfsLowerBound();
				BigInteger unprunedUpper = dee.pmat.calcUnprunedConfsUpperBound();
				BigInteger prunedLower = allConfs.subtract(unprunedUpper);
				BigInteger prunedUpper = allConfs.subtract(unprunedLower);
				double percentPrunedLower = 100.0*prunedLower.doubleValue()/allConfs.doubleValue();
				double percentPrunedUpper = 100.0*prunedUpper.doubleValue()/allConfs.doubleValue();
				log("Conformations defined by conformation space:                           %s", formatBig(allConfs));
				log("Conformations pruned (by singles and pairs, bounds):                   [%s,%s]", formatBig(prunedLower), formatBig(prunedUpper));
				log("Conformations remaining after pruning (by singles and pairs, bounds):  [%s,%s]", formatBig(unprunedLower), formatBig(unprunedUpper));
				log("Percent conformations pruned (by singles and pairs, bounds):           [%.6f,%.6f]", percentPrunedLower, percentPrunedUpper);
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
		pmat.forEachUnprunedSingle((candidatePos, candidateRc) -> {

			// is there a competitor rc that has much lower energy?
			PruningMatrix.IteratorCommand result = pmat.forEachUnprunedSingleAt(candidatePos, (competitorPos, competitorRc) -> {

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
					for (int witnessRc=0; witnessRc<confSpace.positions.get(witnessPos).resConfs.size(); witnessRc++) {

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
		pmat.forEachUnprunedPair((candidatePos1, candidateRc1, candidatePos2, candidateRc2) -> {

			// prune the candidate rc if we can find a competitor rc that has much lower energy
			PruningMatrix.IteratorCommand result = pmat.forEachUnprunedPairAt(candidatePos1, candidatePos2, (competitorPos1, competitorRc1, competitorPos2, competitorRc2) -> {

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
					for (int witnessRc=0; witnessRc<confSpace.positions.get(witnessPos).resConfs.size(); witnessRc++) {

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
				}

				if (energyDiffSum > energyDiffThreshold) {
					return PruningMatrix.IteratorCommand.Break;
				} else {
					return PruningMatrix.IteratorCommand.Continue;
				}
			});

			// if the iteration terminated early (ie, we found a suitable competitor), then prune the candidate
			if (result == PruningMatrix.IteratorCommand.Break) {
				pmat.prunePair(candidatePos1, candidateRc1, candidatePos2, candidateRc2);
			}

			// always go onto to check the next candidate
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	public void pruneTriplesGoldstein(double energyDiffThreshold, boolean typeDependent) {
		pmat.forEachUnprunedTriple((candidatePos1, candidateRc1, candidatePos2, candidateRc2, candidatePos3, candidateRc3) -> {

			// prune the candidate rc if we can find a competitor rc that has much lower energy
			PruningMatrix.IteratorCommand result = pmat.forEachUnprunedTripleAt(candidatePos1, candidatePos2, candidatePos3, (competitorPos1, competitorRc1, competitorPos2, competitorRc2, competitorPos3, competitorRc3) -> {

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
					for (int witnessRc=0; witnessRc<confSpace.positions.get(witnessPos).resConfs.size(); witnessRc++) {

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
				}

				if (energyDiffSum > energyDiffThreshold) {
					return PruningMatrix.IteratorCommand.Break;
				} else {
					return PruningMatrix.IteratorCommand.Continue;
				}
			});

			// if the iteration terminated early (ie, we found a suitable competitor), then prune the candidate
			if (result == PruningMatrix.IteratorCommand.Break) {
				pmat.pruneTriple(candidatePos1, candidateRc1, candidatePos2, candidateRc2, candidatePos3, candidateRc3);
			}

			// always go onto to check the next candidate
			return PruningMatrix.IteratorCommand.Continue;
		});
	}
}
