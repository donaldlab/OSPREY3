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

package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.function.Function;

public class KStarScore {

	public final PartitionFunction.Result protein;
	public final PartitionFunction.Result ligand;
	public final PartitionFunction.Result complex;

	/** the K* score, guaranteed to be an epsilon-approximation or null */
	public final BigDecimal score;

	/** lower bound on K* score */
	public final BigDecimal lowerBound;

	/** upper bound on K* score, could be singleton instance MathTools.BigPositiveInfinity */
	public final BigDecimal upperBound;

	public KStarScore(BigDecimal score, BigDecimal lowerBound, BigDecimal upperBound) {

		this.protein = null;
		this.ligand = null;
		this.complex = null;

		this.score = score;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
	}

	public KStarScore(PartitionFunction.Result protein, PartitionFunction.Result ligand, PartitionFunction.Result complex) {

		this.protein = protein;
		this.ligand = ligand;
		this.complex = complex;

		// calculate the K* score
		// or don't. I'm not the boss of you
		if (protein.status == PartitionFunction.Status.Estimated
				&& ligand.status == PartitionFunction.Status.Estimated
				&& complex.status == PartitionFunction.Status.Estimated) {

			BigDecimal x = MathTools.bigDivideDivide(
				complex.values.qstar,
				protein.values.qstar,
				ligand.values.qstar,
				PartitionFunction.decimalPrecision
			);
			if (MathTools.isNaN(x)) {
				this.score = null;
			} else {
				this.score = x;
			}

		} else {

			// pfuncs not estimated enough, so we don't have a score at all
			this.score = null;
		}

		// calc the lower bound
		this.lowerBound = MathTools.bigDivideDivide(
			complex.values.calcLowerBound(),
			protein.values.calcUpperBound(),
			ligand.values.calcUpperBound(),
			PartitionFunction.decimalPrecision
		);

		// calc the upper bound
		this.upperBound = MathTools.bigDivideDivide(
			complex.values.calcUpperBound(),
			protein.values.calcLowerBound(),
			ligand.values.calcLowerBound(),
			PartitionFunction.decimalPrecision
		);
	}

	public static boolean isLigandComplexUseful(PartitionFunction.Result protein) {

		// assuming we compute pfuncs in order of: protein, ligand, complex:
		// unbound stability is the only thing that would
		// make the ligand or complex pfunc results useless at this point
		return protein.status != PartitionFunction.Status.Unstable;
	}

	public static boolean isComplexUseful(PartitionFunction.Result protein, PartitionFunction.Result ligand) {

		// assuming we compute pfuncs in order of: protein, ligand, complex:
		// unbound stability is the only thing that would
		// make the complex pfunc results useless at this point
		return protein.status != PartitionFunction.Status.Unstable
			&& ligand.status != PartitionFunction.Status.Unstable;
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
		return String.format("%-9s in [%-9s,%9s] (log10)",
			trim.apply(scoreLog10String()),
			trim.apply(lowerBoundLog10String()),
			trim.apply(upperBoundLog10String())
		);
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof KStarScore && equals((KStarScore)other);
	}

	public boolean equals(KStarScore other) {
		return MathTools.isSameValue(this.score, other.score)
			&& MathTools.isSameValue(this.lowerBound, other.lowerBound)
			&& MathTools.isSameValue(this.upperBound, other.upperBound);
	}

	public boolean isSimilarTo(KStarScore other, double relativeEpsilon) {
		return MathTools.isRelativelySame(this.score, other.score, PartitionFunction.decimalPrecision, relativeEpsilon)
			&& MathTools.isRelativelySame(this.lowerBound, other.lowerBound, PartitionFunction.decimalPrecision, relativeEpsilon)
			&& MathTools.isRelativelySame(this.upperBound, other.upperBound, PartitionFunction.decimalPrecision, relativeEpsilon);
	}

	public Double scoreLog10() { return scoreToLog10(score); }
	public Double lowerBoundLog10() { return scoreToLog10(lowerBound); }
	public Double upperBoundLog10() { return scoreToLog10(upperBound); }

	public String scoreLog10String() { return scoreToLog10String(score); }
	public String lowerBoundLog10String() { return scoreToLog10String(lowerBound); }
	public String upperBoundLog10String() { return scoreToLog10String(upperBound); }

	public static Double scoreToLog10(BigDecimal score) {
		if (score != null) {
			return Math.log10(score.doubleValue());
		}
		return null;
	}

	public static String scoreToString(BigDecimal score) {
		if (score != null) {
			return String.format("%e", score.doubleValue());
		}
		return "none";
	}

	public static String scoreToLog10String(BigDecimal score) {
		if (score != null) {
			return String.format("%f", scoreToLog10(score));
		}
		return "none";
	}
}
