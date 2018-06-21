package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.function.Function;

public class EWAKStarScore {

	public final EWAKStarPartitionFunction.Result protein;
	public final EWAKStarPartitionFunction.Result ligand;
	public final EWAKStarPartitionFunction.Result complex;

	/** the K* score, guaranteed to be an epsilon-approximation or null */
	public final BigDecimal score;

	/** lower bound on K* score */
	public final BigDecimal lowerBound;

	/** upper bound on K* score, could be singleton instance MathTools.BigPositiveInfinity */
	public final BigDecimal upperBound;

	public EWAKStarScore(BigDecimal score, BigDecimal lowerBound, BigDecimal upperBound) {

		this.protein = null;
		this.ligand = null;
		this.complex = null;

		this.score = score;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
	}

	public EWAKStarScore(EWAKStarGradientDescentPfunc.Result protein, EWAKStarGradientDescentPfunc.Result ligand, EWAKStarGradientDescentPfunc.Result complex) {

		this.protein = protein;
		this.ligand = ligand;
		this.complex = complex;

		// calculate the K* score
		// or don't. I'm not the boss of you
		if (protein.status == EWAKStarPartitionFunction.Status.ConfLimitReached
				&& ligand.status == EWAKStarPartitionFunction.Status.ConfLimitReached
				&& complex.status == EWAKStarPartitionFunction.Status.ConfLimitReached ||
				protein.status == EWAKStarPartitionFunction.Status.EpsilonReached
						&& ligand.status == EWAKStarPartitionFunction.Status.EpsilonReached
						&& complex.status == EWAKStarPartitionFunction.Status.EpsilonReached ||
				protein.status == EWAKStarPartitionFunction.Status.EnergyReached
						&& ligand.status == EWAKStarPartitionFunction.Status.EnergyReached
						&& complex.status == EWAKStarPartitionFunction.Status.EnergyReached) {

			BigDecimal x = MathTools.bigDivideDivide(
				complex.values.qstar,
				protein.values.qstar,
				ligand.values.qstar,
				EWAKStarGradientDescentPfunc.decimalPrecision
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
			EWAKStarGradientDescentPfunc.decimalPrecision
		);

		// calc the upper bound
		this.upperBound = MathTools.bigDivideDivide(
			complex.values.calcUpperBound(),
			protein.values.calcLowerBound(),
			ligand.values.calcLowerBound(),
			EWAKStarGradientDescentPfunc.decimalPrecision
		);
	}

	public static boolean isLigandComplexUseful(EWAKStarGradientDescentPfunc.Result protein) {

		// assuming we compute pfuncs in order of: protein, ligand, complex:
		// unbound stability is the only thing that would
		// make the ligand or complex pfunc results useless at this point
		return protein.status != EWAKStarGradientDescentPfunc.Status.Unstable;
	}

	public static boolean isComplexUseful(EWAKStarGradientDescentPfunc.Result protein, EWAKStarGradientDescentPfunc.Result ligand) {

		// assuming we compute pfuncs in order of: protein, ligand, complex:
		// unbound stability is the only thing that would
		// make the complex pfunc results useless at this point
		return protein.status != EWAKStarGradientDescentPfunc.Status.Unstable
			&& ligand.status != EWAKStarGradientDescentPfunc.Status.Unstable;
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
		return other instanceof EWAKStarScore && equals((EWAKStarScore)other);
	}

	public boolean equals(EWAKStarScore other) {
		return MathTools.isSameValue(this.score, other.score)
			&& MathTools.isSameValue(this.lowerBound, other.lowerBound)
			&& MathTools.isSameValue(this.upperBound, other.upperBound);
	}

	public boolean isSimilarTo(EWAKStarScore other, double relativeEpsilon) {
		return MathTools.isRelativelySame(this.score, other.score, EWAKStarGradientDescentPfunc.decimalPrecision, relativeEpsilon)
			&& MathTools.isRelativelySame(this.lowerBound, other.lowerBound, EWAKStarGradientDescentPfunc.decimalPrecision, relativeEpsilon)
			&& MathTools.isRelativelySame(this.upperBound, other.upperBound, EWAKStarGradientDescentPfunc.decimalPrecision, relativeEpsilon);
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
