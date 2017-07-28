package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.math.BigDecimal;
import java.math.RoundingMode;

public class KStarScore {

	public final PartitionFunction.Result protein;
	public final PartitionFunction.Result ligand;
	public final PartitionFunction.Result complex;

	/** guaranteed to be an epsilon-approximation to K* or null */
	public final BigDecimal score;

	public final BigDecimal lowerBound;
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

		// calculate the lower bound
		BigDecimal proteinUpperBound = protein.values.calcUpperBound();
		BigDecimal ligandUpperBound = ligand.values.calcUpperBound();
		if (proteinUpperBound.compareTo(BigDecimal.ZERO) == 0 || ligandUpperBound.compareTo(BigDecimal.ZERO) == 0) {
			// can't represent infinity in BigDecimal, so use null instead
			this.lowerBound = null;
		} else {
			this.lowerBound = complex.values.calcLowerBound()
				.divide(proteinUpperBound, RoundingMode.HALF_UP)
				.divide(ligandUpperBound, RoundingMode.HALF_UP);
		}

		// calculate the upper bound
		BigDecimal proteinLowerBound = protein.values.calcLowerBound();
		BigDecimal ligandLowerBound = ligand.values.calcLowerBound();
		if (proteinLowerBound.compareTo(BigDecimal.ZERO) == 0 || ligandLowerBound.compareTo(BigDecimal.ZERO) == 0) {
			// can't represent infinity in BigDecimal, so use null instead
			this.upperBound = null;
		} else {
			this.upperBound = complex.values.calcUpperBound()
				.divide(proteinLowerBound, RoundingMode.HALF_UP)
				.divide(ligandLowerBound, RoundingMode.HALF_UP);
		}
	}

	@Override
	public String toString() {
		return String.format("%s in [%s,%s]",
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
