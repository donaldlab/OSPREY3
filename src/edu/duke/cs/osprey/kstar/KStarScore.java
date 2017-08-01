package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.math.RoundingMode;

import static edu.duke.cs.osprey.tools.MathTools.isZero;

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

	public KStarScore(PartitionFunction.Result protein, PartitionFunction.Result ligand, PartitionFunction.Result complex) {

		this.protein = protein;
		this.ligand = ligand;
		this.complex = complex;

		// calculate the K* score
		if (protein.status == PartitionFunction.Status.Estimated
				&& ligand.status == PartitionFunction.Status.Estimated
				&& complex.status == PartitionFunction.Status.Estimated) {

			this.score = complex.values.qstar
				.divide(protein.values.qstar, RoundingMode.HALF_UP)
				.divide(ligand.values.qstar, RoundingMode.HALF_UP);
		} else {
			this.score = null;
		}

		// calculate the K* bounds
		BigDecimal proteinLowerBound = protein.values.calcLowerBound();
		BigDecimal proteinUpperBound = protein.values.calcUpperBound();
		BigDecimal ligandLowerBound = ligand.values.calcLowerBound();
		BigDecimal ligandUpperBound = ligand.values.calcUpperBound();
		BigDecimal complexLowerBound = complex.values.calcLowerBound();
		BigDecimal complexUpperBound = complex.values.calcUpperBound();

		if (isZero(proteinUpperBound) || isZero(ligandUpperBound) || isZero(complexUpperBound)) {

			// these must be highly unfavorable structures
			// Boltzmann says they will happen only with infinitesimally small probability
			// but we apparently ran out of precision calculating e^(-energy), so let's just call it zero
			this.lowerBound = BigDecimal.ZERO;
			this.upperBound = BigDecimal.ZERO;

		} else {

			this.lowerBound = complexLowerBound
				.divide(proteinUpperBound, RoundingMode.HALF_UP)
				.divide(ligandUpperBound, RoundingMode.HALF_UP);

			if (isZero(proteinLowerBound) || isZero(ligandLowerBound)) {
				// this *should* be impossible for a single-sequence bound (since if the lower bound is zero,
				// the upper bound *should* be zero too), but it could easily happen for a multi-sequence bound
				this.upperBound = MathTools.BigPositiveInfinity;
			} else {
				this.upperBound = complexUpperBound
					.divide(proteinLowerBound, RoundingMode.HALF_UP)
					.divide(ligandLowerBound, RoundingMode.HALF_UP);
			}
		}
	}

	@Override
	public String toString() {
		return String.format("%-9s in [%-9s,%-9s]",
			scoreLog10String(),
			lowerBoundLog10String(),
			upperBoundLog10String()
		);
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
