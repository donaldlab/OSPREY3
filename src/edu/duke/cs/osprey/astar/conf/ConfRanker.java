package edu.duke.cs.osprey.astar.conf;

import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;

import java.math.BigInteger;
import java.util.Arrays;

public class ConfRanker {

	public static interface ConfAStarConfigurer {
		void configure(ConfAStarTree.Builder builder);
	}

	private static final int NotAssigned = -1;

	public final SimpleConfSpace confSpace;
	public final EnergyMatrix emat;
	public final ConfAStarConfigurer astarConfigurer;

	private final EnergyMatrix negatedEmat;
	private final int[] noAssignmentsMask;
	private final AStarScorer scorer;
	private final RCs allRCs;
	private final ConfIndex confIndex;

	public ConfRanker(SimpleConfSpace confSpace, EnergyMatrix emat, ConfAStarConfigurer astarConfigurer) {

		this.confSpace = confSpace;
		this.emat = emat;
		this.astarConfigurer = astarConfigurer;

		negatedEmat = new NegatedEnergyMatrix(confSpace, emat);

		// make the all-unassigned conf mask
		noAssignmentsMask = new int[confSpace.positions.size()];
		Arrays.fill(noAssignmentsMask, NotAssigned);

		// get all the RCs in the conf space
		allRCs = new RCs(confSpace);

		// get the A* scorer
		scorer = makeAstar(emat, allRCs).gscorer;

		// allocate the conf index for conf scoring, and prep for fully-assigned confs
		confIndex = new ConfIndex(confSpace.positions.size());
		confIndex.numDefined = confSpace.positions.size();
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			 confIndex.definedPos[pos.index] = pos.index;
		}
		confIndex.numUndefined = 0;
	}

	private ConfAStarTree makeAstar(EnergyMatrix emat, RCs rcs) {
		ConfAStarTree.Builder builder = new ConfAStarTree.Builder(emat, rcs);
		astarConfigurer.configure(builder);
		return builder.build();
	}

	public BigInteger getNumConfsAtMost(double queryScore) {
		return getNumConfsAtMost(queryScore, noAssignmentsMask);
	}

	private BigInteger getNumConfsAtMost(double queryScore, int[] confMask) {

		// NOTE: conf masks are always sequences of assigned positions, then unassigned positions
		// e.g., -1, -1, -1
		// or 0, 5, -1
		// but never -1, -1, 2

		// find the first unassigned position, if any
		SimpleConfSpace.Position firstUnassignedPos = confSpace.positions.stream()
			.filter((pos) -> confMask[pos.index] == -1)
			.findFirst()
			.orElse(null);

		// if all positions are assigned (and confMask is really just a conf),
		// then just compare conf scores directly
		if (firstUnassignedPos == null) {

			// update the conf index with the conf RCs
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				confIndex.definedRCs[pos.index] = confMask[pos.index];
			}
			double confScore = scorer.calc(confIndex, allRCs);

			if (confScore <= queryScore) {
				return BigInteger.ONE;
			} else {
				return BigInteger.ZERO;
			}
		}

		// get the RCs for this sub-tree using the conf mask
		RCs rcs = new RCs(confSpace, (pos, resConf) -> {
			int rc = confMask[pos.index];
			return rc == -1 || rc == resConf.index;
		});

		// no flexibility? no conformations
		if (!rcs.hasConfs()) {
			return BigInteger.ZERO;
		}

		// calculate the min and max scores on this sub-tree
		// TODO: could do this in parallel?
		double minScore = makeAstar(emat, rcs)
			.nextConf()
			.getScore();
		double maxScore = -makeAstar(negatedEmat, rcs)
			.nextConf()
			.getScore();

		// compare the sub-tree bounds to the query score
		if (minScore > queryScore) {

			// the whole sub-tree is above the score,
			// so no confs are at most the score
			return BigInteger.ZERO;

		} else if (maxScore <= queryScore) {

			// the whole sub-tree is at most the score
			return rcs.getNumConformations();
		}

		// the sub-tree is part below and part above, so recurse on the sub-sub-trees
		// TODO: could do this in parallel?
		BigInteger numConfs = BigInteger.ZERO;
		int [] subConfMask = confMask.clone();
		for (SimpleConfSpace.ResidueConf resConf : firstUnassignedPos.resConfs) {

			subConfMask[firstUnassignedPos.index] = resConf.index;

			// recurse
			BigInteger subNumConfs = getNumConfsAtMost(queryScore, subConfMask);

			numConfs = numConfs.add(subNumConfs);
		}

		return numConfs;
	}
}
