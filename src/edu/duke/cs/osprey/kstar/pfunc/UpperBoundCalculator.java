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

package edu.duke.cs.osprey.kstar.pfunc;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.math.BigInteger;

public class UpperBoundCalculator {

	public final ConfSearch tree;
	public BigInteger numUnscoredConfs;

	public int numScoredConfs = 0;
	public BigDecimal weightedScoreSum = BigDecimal.ZERO;
	public BigDecimal unscoredBound = BigDecimal.ZERO;
	public BigDecimal totalBound = BigDecimal.ZERO;
	public double delta = 1.0;

	private BoltzmannCalculator boltzmann = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

	public UpperBoundCalculator(ConfSearch tree, BigInteger numConfsBeforePruning) {
		this.tree = tree;
		this.numUnscoredConfs = numConfsBeforePruning;
	}

	public UpperBoundCalculator run(int numConfs) {
		boolean canContinue = true;
		while (canContinue && numConfs > 0) {
			canContinue = scoreNextConf();
			numConfs--;
		}
		return this;
	}

	public UpperBoundCalculator run(int numConfs, double epsilon) {
		boolean canContinue = delta > epsilon;
		while (canContinue && numConfs > 0) {
			canContinue = scoreNextConf() && delta > epsilon;
			numConfs--;
		}
		return this;
	}

	public boolean scoreNextConf() {

		// get the next conf
		ConfSearch.ScoredConf conf = tree.nextConf();
		if (conf == null) {
			return false;
		}

		numScoredConfs++;

		// compute the boltzmann weight for this conf
		BigDecimal weightedScore = boltzmann.calc(conf.getScore());

		// update counters/sums/bounds
		weightedScoreSum = weightedScoreSum.add(weightedScore);
		numUnscoredConfs = numUnscoredConfs.subtract(BigInteger.ONE);
		unscoredBound = weightedScore.multiply(new BigDecimal(numUnscoredConfs));
		totalBound = weightedScoreSum.add(unscoredBound);

		// update delta
		delta = MathTools.bigDivide(unscoredBound, totalBound, PartitionFunction.decimalPrecision).doubleValue();

		// keep going if the boltzmann weight is greater than zero
		return MathTools.isGreaterThan(weightedScore, BigDecimal.ZERO);
	}

	@Override
	public String toString() {
		return String.format("UpperBoundCalculator  scored: %8d  sum: %e   unscored: %e  %e   bound: %e   delta: %f",
			numScoredConfs,
			weightedScoreSum.doubleValue(),
			numUnscoredConfs.doubleValue(),
			unscoredBound.doubleValue(),
			totalBound.doubleValue(),
			delta
		);
	}
}
