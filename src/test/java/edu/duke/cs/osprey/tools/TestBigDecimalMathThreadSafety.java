package edu.duke.cs.osprey.tools;

import static edu.duke.cs.osprey.TestBase.isRelatively;
import static org.hamcrest.MatcherAssert.*;

import ch.obermuhlner.math.big.BigDecimalMath;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import org.junit.jupiter.api.Test;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;


/**
 * These tests won't always fail if there's a problem!
 * Concurrency bugs are awful like that. T_T
 * But if they do fail, you know for sure there's a problem.
 */
public class TestBigDecimalMathThreadSafety {

	@Test
	public void exp() {

		var es = new double[] {
			-1524.0681,
			-1523.9364,
			-1524.2018,
			-1523.7533
		};

		var zs = new BigDecimal[es.length];

		var mc = new MathContext(8, RoundingMode.HALF_UP);

		// compute all the exp values in parallel
		try (var tasks = new ThreadPoolTaskExecutor()) {
			tasks.start(es.length);

			for (int i=0; i<es.length; i++) {
				final int fi = i;
				tasks.submit(
					() -> BigDecimalMath.exp(BigDecimal.valueOf(-es[fi]), mc),
					z -> zs[fi] = z
				);
			}
			tasks.waitForFinish();
		}

		// check the results
		final double epsilon = 1e-8;
		assertThat(zs[0], isRelatively(new BigDecimal("7.8408993E+661"), epsilon));
		assertThat(zs[1], isRelatively(new BigDecimal("6.8733632E+661"), epsilon));
		assertThat(zs[2], isRelatively(new BigDecimal("8.9625388E+661"), epsilon));
		assertThat(zs[3], isRelatively(new BigDecimal("5.7233456E+661"), epsilon));
	}
}
