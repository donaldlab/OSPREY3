package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;
import org.hamcrest.BaseMatcher;
import org.hamcrest.Description;
import org.hamcrest.Matcher;
import org.junit.Test;

import java.math.BigDecimal;

public class TestKStarScore {

	public static PartitionFunction.Result makeResult(PartitionFunction.Status status, double min, double max) {
		PartitionFunction.Values values = new PartitionFunction.Values();
		values.qstar = BigDecimal.valueOf(min);
		if (max == Double.POSITIVE_INFINITY) {
			values.qprime = MathTools.BigPositiveInfinity;
		} else {
			values.qprime = BigDecimal.valueOf(max - min);
		}
		return new PartitionFunction.Result(status, values, 0);
	}

	public static KStarScore makeScore(Double score, double min, double max) {
		return new KStarScore(
			score != null ? MathTools.biggen(score) : null,
			MathTools.biggen(min),
			MathTools.biggen(max)
		);
	}

	public static Matcher<KStarScore> isRelatively(final KStarScore expected) {
		return new BaseMatcher<KStarScore>() {

			@Override
			public boolean matches(Object obj) {
				final double epsilon = 1e-8;
				return expected.isSimilarTo((KStarScore)obj, epsilon);
			}

			@Override
			public void describeTo(Description desc) {
				desc.appendText("close to ").appendValue(expected);
			}

			@Override
			public void describeMismatch(Object obj, Description desc) {
				desc.appendText("was      ").appendValue(obj);
			}
		};
	}

	@Test
	public void resultFactory() {
		PartitionFunction.Result result = makeResult(PartitionFunction.Status.Estimated, 5, 9);
		assertThat(result.values.qstar.doubleValue(), is(5.0));
		assertThat(result.values.qprime.doubleValue(), is(4.0));
		assertThat(result.values.calcLowerBound().doubleValue(), is(5.0));
		assertThat(result.values.calcUpperBound().doubleValue(), is(9.0));
	}


	// stability tests

	@Test // PLC
	public void allStable() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(4.0/3/2, 4.0/7/6, 8.0/3/2)));
	}

	@Test // pLC
	public void onlyProteinUnstable() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)));
	}

	@Test // PlC
	public void onlyLigandUnstable() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)));
	}

	@Test // PLc
	public void onlyComplexUnstable() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0)
		), isRelatively(makeScore(0.0, 0.0, 0.0)));
	}

	@Test // plC
	public void proteinAndLigandUnstable() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)));
	}

	@Test // pLc
	public void proteinAndComplexUnstable() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0)
		), isRelatively(makeScore(null, Double.NaN, Double.NaN)));
	}

	@Test // Plc
	public void ligandAndComplexUnstable() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0)
		), isRelatively(makeScore(null, Double.NaN, Double.NaN)));
	}

	@Test // plc
	public void allUnstable() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 0.0)
		), isRelatively(makeScore(null, Double.NaN, Double.NaN)));
	}


	// lower bound tests

	@Test // pLC
	public void proteinZeroLower() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(Double.POSITIVE_INFINITY, 4.0/7/6, Double.POSITIVE_INFINITY)));
	}

	@Test // PlC
	public void ligandZeroLower() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(Double.POSITIVE_INFINITY, 4.0/7/6, Double.POSITIVE_INFINITY)));
	}

	@Test // PLc
	public void complexZeroLower() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 8.0)
		), isRelatively(makeScore(0.0, 0.0, 8.0/3/2)));
	}

	@Test // plC
	public void proteinAndLigandZeroLower() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(Double.POSITIVE_INFINITY, 4.0/7/6, Double.POSITIVE_INFINITY)));
	}

	@Test // pLc
	public void proteinAndComplexZeroLower() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 8.0)
		), isRelatively(makeScore(null, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // Plc
	public void ligandAndComplexZeroLower() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 8.0)
		), isRelatively(makeScore(null, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // plc
	public void allZeroLower() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, 8.0)
		), isRelatively(makeScore(null, 0.0, Double.POSITIVE_INFINITY)));
	}


	// upper bound tests

	@Test // pLC
	public void proteinInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(4.0/3/2, 0.0, 8.0/3/2)));
	}

	@Test // PlC
	public void ligandInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(4.0/3/2, 0.0, 8.0/3/2)));
	}

	@Test // PLc
	public void complexInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, Double.POSITIVE_INFINITY)
		), isRelatively(makeScore(4.0/3/2, 4.0/7/6, Double.POSITIVE_INFINITY)));
	}

	@Test // plC
	public void proteinAndLigandInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 3.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(4.0/3/2, 0.0, 8.0/3/2)));
	}

	@Test // pLc
	public void proteinAndComplexInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, Double.POSITIVE_INFINITY)
		), isRelatively(makeScore(4.0/3/2, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // Plc
	public void ligandAndComplexInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 4.0, Double.POSITIVE_INFINITY)
		), isRelatively(makeScore(4.0/3/2, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // plc
	public void allInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 3.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 4.0, Double.POSITIVE_INFINITY)
		), isRelatively(makeScore(4.0/3/2, 0.0, Double.POSITIVE_INFINITY)));
	}


	// lower and upper bound tests

	@Test // pLC
	public void proteinZeroLowerInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(Double.POSITIVE_INFINITY, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // PlC
	public void ligandZeroLowerInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(Double.POSITIVE_INFINITY, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // PLc
	public void complexZeroLowerInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY)
		), isRelatively(makeScore(0.0, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // plC
	public void proteinAndLigandZeroLowerInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(Double.POSITIVE_INFINITY, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // pLc
	public void proteinAndComplexZeroLowerInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY)
		), isRelatively(makeScore(null, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // Plc
	public void ligandAndComplexZeroLowerInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY)
		), isRelatively(makeScore(null, 0.0, Double.POSITIVE_INFINITY)));
	}

	@Test // plc
	public void allZeroLowerInfUpper() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY),
			makeResult(PartitionFunction.Status.Estimated, 0.0, Double.POSITIVE_INFINITY)
		), isRelatively(makeScore(null, 0.0, Double.POSITIVE_INFINITY)));
	}


	// not estimated tests

	@Test
	public void proteinNotEstimated() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimating, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(null, 4.0/7/6, 8.0/3/2)));
	}

	@Test
	public void ligandNotEstimated() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimating, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimated, 4.0, 8.0)
		), isRelatively(makeScore(null, 4.0/7/6, 8.0/3/2)));
	}

	@Test
	public void complexNotEstimated() {
		assertThat(new KStarScore(
			makeResult(PartitionFunction.Status.Estimated, 2.0, 6.0),
			makeResult(PartitionFunction.Status.Estimated, 3.0, 7.0),
			makeResult(PartitionFunction.Status.Estimating, 4.0, 8.0)
		), isRelatively(makeScore(null, 4.0/7/6, 8.0/3/2)));
	}
}
