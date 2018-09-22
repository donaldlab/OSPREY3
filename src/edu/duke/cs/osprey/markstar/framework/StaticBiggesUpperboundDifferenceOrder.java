package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.*;

public class StaticBiggesUpperboundDifferenceOrder implements edu.duke.cs.osprey.astar.conf.order.AStarOrder {

	public final MathTools.Optimizer optimizer;
	private List<Integer> posOrder;
	private BoltzmannCalculator calculator = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

	public StaticBiggesUpperboundDifferenceOrder() {
		this(MathTools.Optimizer.Maximize);
	}

	public StaticBiggesUpperboundDifferenceOrder(MathTools.Optimizer optimizer) {
		this.optimizer = optimizer;
	}

	private AStarScorer gscorer;
	private AStarScorer hscorer;

	@Override
	public void setScorers(AStarScorer gscorer, AStarScorer hscorer) {
		this.gscorer = gscorer;
		this.hscorer = hscorer;
	}

	@Override
	public boolean isDynamic() {
		return false;
	}

	@Override
	public int getNextPos(ConfIndex confIndex, RCs rcs) {
		if (posOrder == null) {
			posOrder = calcPosOrder(confIndex, rcs);
		}
		return posOrder.get(confIndex.node.getLevel());
	}

	private List<Integer> calcPosOrder(ConfIndex confIndex, RCs rcs) {
		// init permutation array with only undefined positions and score them
		List<Integer> undefinedOrder = new ArrayList<Integer>();
		Map<Integer,BigDecimal> scores = new TreeMap<>();
		for (int posi=0; posi<confIndex.numUndefined; posi++) {
			int pos = confIndex.undefinedPos[posi];
			undefinedOrder.add(pos);
			scores.put(pos, scorePos(confIndex, rcs, pos));
		}

		// sort positions in order of decreasing score
		Collections.sort(undefinedOrder, new Comparator<Integer>() {

			@Override
			public int compare(Integer pos1, Integer pos2) {
				BigDecimal score1 = scores.get(pos1);
				BigDecimal score2 = scores.get(pos2);
				// NOTE: use reverse order for decreasing sort
				return score2.compareTo(score1);
			}
		});

		// prepend the defined positions to build the full order
		List<Integer> order = new ArrayList<>();
		for (int posi=0; posi<confIndex.numDefined; posi++) {
			int pos = confIndex.definedPos[posi];
			order.add(pos);
		}
		order.addAll(undefinedOrder);


		return order;
	}

	private void putResidueInPosition(List<Integer> order, int residueIndex, int orderIndex) {
	    swapPos(order, residueIndex, order.get(orderIndex));
	}

	private void swapPos(List<Integer> order, int a, int b) {
		int aIndex = order.indexOf(a);
		int bIndex = order.indexOf(b);

		order.set(bIndex, a);
		order.set(aIndex, b);
	}


	BigDecimal scorePos(ConfIndex confIndex, RCs rcs, int pos) {

		// check all the RCs at this pos and aggregate the energies
		double parentScore = confIndex.node.getScore();
		double reciprocalSum = 0;
		BigDecimal maxUpper = MathTools.BigNegativeInfinity;
		BigDecimal minUpper = MathTools.BigPositiveInfinity;
		for (int rc : rcs.get(pos)) {
			double childScore = gscorer.calcDifferential(confIndex, rcs, pos, rc)
				+ hscorer.calcDifferential(confIndex, rcs, pos, rc);
			BigDecimal childWeightedScore =  calculator.calc(childScore);
			if(MathTools.isLessThan(childWeightedScore, minUpper))
				minUpper = childWeightedScore;
			if(MathTools.isGreaterThan(childWeightedScore, maxUpper))
				maxUpper = childWeightedScore;
		}

		return maxUpper.subtract(minUpper);
	}
}
