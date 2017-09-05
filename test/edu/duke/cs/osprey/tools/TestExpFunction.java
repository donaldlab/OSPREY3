package edu.duke.cs.osprey.tools;

import org.junit.Test;

import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.*;

public class TestExpFunction {

	@Test
	public void exp() {
		final double Epsilon = 1e-8;
		ExpFunction exp = new ExpFunction();
		assertThat(exp.exp(-10.0).doubleValue(), isAbsolutely(Math.exp(-10.0), Epsilon));
		assertThat(exp.exp( -5.0).doubleValue(), isAbsolutely(Math.exp( -5.0), Epsilon));
		assertThat(exp.exp( -1.0).doubleValue(), isAbsolutely(Math.exp( -1.0), Epsilon));
		assertThat(exp.exp(  0.0).doubleValue(), isAbsolutely(Math.exp(  0.0), Epsilon));
		assertThat(exp.exp(  1.0).doubleValue(), isAbsolutely(Math.exp(  1.0), Epsilon));
		assertThat(exp.exp(  5.0).doubleValue(), isAbsolutely(Math.exp(  5.0), Epsilon));
		assertThat(exp.exp( 10.0).doubleValue(), isAbsolutely(Math.exp( 10.0), Epsilon));
	}
}