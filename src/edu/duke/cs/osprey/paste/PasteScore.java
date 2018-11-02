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

package edu.duke.cs.osprey.paste;

import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.function.Function;

public class PasteScore {

	public static double constRT = -1.9891/1000.0 * 298.15;//RT in kcal/mol
	public final PastePartitionFunction.Result protein;
	public final PastePartitionFunction.Result wt;

	/** the K* score, guaranteed to be an epsilon-approximation or null */
	public final Double score;

	/** lower bound on K* score */
	public final Double lowerBound;

	/** upper bound on K* score, could be singleton instance MathTools.BigPositiveInfinity */
	public final Double upperBound;

	public String stability;

	public PasteScore(Double score, Double lowerBound, Double upperBound) {

		this.protein = null;
		this.wt = null;

		this.score = score;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
	}

	public PasteScore(PastePartitionFunction.Result protein, PastePartitionFunction.Result wt) {

		this.protein = protein;
		this.wt = wt;

		// report the partition function value
		if (protein.status == PastePartitionFunction.Status.EnergyReached ||
				protein.status == PastePartitionFunction.Status.ConfLimitReached ||
				protein.status == PastePartitionFunction.Status.EpsilonReached ||
				protein.status == PastePartitionFunction.Status.NoWindowOverlap) {

			if (MathTools.isNaN(protein.values.qstar)) {
				this.score = null;
			} else {
				this.score = constRT*Math.log(new BigMath(PastePartitionFunction.decimalPrecision)
						.set(protein.values.calcLowerBound())
						.div(wt.values.calcLowerBound())
						.get().doubleValue());
			}

		} else {

			// pfuncs not estimated enough, so we don't have a score at all
			this.score = null;
		}

		// calc the lower bound
		this.lowerBound = constRT*Math.log(new BigMath(PastePartitionFunction.decimalPrecision)
				.set(protein.values.calcUpperBound())
				.div(wt.values.calcLowerBound())
				.get().doubleValue());

		// calc the upper bound
		this.upperBound = constRT*Math.log((new BigMath(PastePartitionFunction.decimalPrecision)
				.set(protein.values.calcLowerBound())
				.div(wt.values.calcUpperBound())
				.get().doubleValue()));

		if (this.upperBound < 0 && this.lowerBound < 0){
			this.stability = "Mutation Increases Stability";
		} else if (this.upperBound > 0 && this.lowerBound > 0){
			this.stability = "Mutation Decreases Stability";
		} else
			this.stability = "Affect on Stability Unclear";
	}

	@Override
	public String toString() {
		Function<String,String> trim = (s) -> {
			if (s.length() > 9) {
				return s.substring(0, 9);
			} else {
				return s;
			}
		};
		return String.format("%-9s in [%-9s,%9s]       %35s",
			trim.apply(scoreToString(score)),
			trim.apply(scoreToString(lowerBound)),
			trim.apply(scoreToString(upperBound)),
			stability
		);
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof PasteScore && equals((PasteScore)other);
	}

	public boolean equals(PasteScore other) {
		return this.score.equals(other.score)
			&& this.lowerBound.equals(other.lowerBound)
			&& this.upperBound.equals(other.upperBound);
	}

	public String lowerBoundLog10String() { return scoreToLog10String(protein.values.calcLowerBound().doubleValue()); }
	public String upperBoundLog10String() { return scoreToLog10String(protein.values.calcUpperBound().doubleValue()); }

	public static Double scoreToLog10(Double score) {
		if (score != null) {
			return Math.log10(score.doubleValue());
		}
		return null;
	}

	public static String scoreToString(Double score) {
		if (score != null) {
			return String.format("%f", score);
		}
		return "none";
	}

	public static String scoreToLog10String(Double score) {
		if (score != null) {
			return String.format("%f", scoreToLog10(score));
		}
		return "none";
	}
}
