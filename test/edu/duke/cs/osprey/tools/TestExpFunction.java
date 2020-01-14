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

package edu.duke.cs.osprey.tools;

import org.junit.Test;
import org.junit.experimental.theories.suppliers.TestedOn;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.*;

public class TestExpFunction {

	@Test
	public void exp() {
		final double Epsilon = 1e-8;
		ExpFunction exp = new ExpFunction();
		assertThat(exp.exp(-10.0).doubleValue(), isAbsolutely(Math.exp(-10.0), Epsilon));
		assertThat(exp.exp( -7.5).doubleValue(), isAbsolutely(Math.exp( -7.5), Epsilon));
		assertThat(exp.exp( -5.0).doubleValue(), isAbsolutely(Math.exp( -5.0), Epsilon));
		assertThat(exp.exp( -2.5).doubleValue(), isAbsolutely(Math.exp( -2.5), Epsilon));
		assertThat(exp.exp( -1.0).doubleValue(), isAbsolutely(Math.exp( -1.0), Epsilon));
		assertThat(exp.exp( -0.5).doubleValue(), isAbsolutely(Math.exp( -0.5), Epsilon));
		assertThat(exp.exp(  0.0).doubleValue(), isAbsolutely(Math.exp(  0.0), Epsilon));
		assertThat(exp.exp(  0.5).doubleValue(), isAbsolutely(Math.exp(  0.5), Epsilon));
		assertThat(exp.exp(  1.0).doubleValue(), isAbsolutely(Math.exp(  1.0), Epsilon));
		assertThat(exp.exp(  2.5).doubleValue(), isAbsolutely(Math.exp(  2.5), Epsilon));
		assertThat(exp.exp(  5.0).doubleValue(), isAbsolutely(Math.exp(  5.0), Epsilon));
		assertThat(exp.exp(  7.5).doubleValue(), isAbsolutely(Math.exp(  7.5), Epsilon));
		assertThat(exp.exp( 10.0).doubleValue(), isAbsolutely(Math.exp( 10.0), Epsilon));
	}

	@Test
	public void testLog(){
		final double Epsilon = 1e-5;
		final MathContext smallest = new MathContext(8, RoundingMode.HALF_UP);
		final MathContext middle = new MathContext(64, RoundingMode.HALF_UP);
		final MathContext ourDefault = new MathContext(100, RoundingMode.HALF_UP);

		ExpFunction exp = new ExpFunction();
		assertThat(exp.log(exp.exp(-10)).doubleValue(), isAbsolutely(-10.0, Epsilon));
		assertThat(exp.log(exp.exp(-7.5)).doubleValue(), isAbsolutely(-7.5, Epsilon));
		assertThat(exp.log(exp.exp(-5.0)).doubleValue(), isAbsolutely(-5.0, Epsilon));
		assertThat(exp.log(exp.exp(-2.5)).doubleValue(), isAbsolutely(-2.5, Epsilon));
		assertThat(exp.log(exp.exp(-1.0)).doubleValue(), isAbsolutely(-1.0, Epsilon));
		assertThat(exp.log(exp.exp(-0.5)).doubleValue(), isAbsolutely(-0.5, Epsilon));
		assertThat(exp.log(exp.exp(0.0)).doubleValue(), isAbsolutely(0.0, Epsilon));
		assertThat(exp.log(exp.exp(0.5)).doubleValue(), isAbsolutely(0.5, Epsilon));
		assertThat(exp.log(exp.exp(1.0)).doubleValue(), isAbsolutely(1.0, Epsilon));
		assertThat(exp.log(exp.exp(2.5)).doubleValue(), isAbsolutely(2.5, Epsilon));
		assertThat(exp.log(exp.exp(5.0)).doubleValue(), isAbsolutely(5.0, Epsilon));
		assertThat(exp.log(exp.exp(7.5)).doubleValue(), isAbsolutely(7.5, Epsilon));
		assertThat(exp.log(exp.exp(10.0)).doubleValue(), isAbsolutely(10.0, Epsilon));

		// this log function will fail if the value's exponent is much greater than the precision + 36, empirically
		// testing the effect of precision on the logs
		assertNotEquals(exp.log(new BigDecimal(6.091e71, smallest)).doubleValue(), isAbsolutely(165.290353874, Epsilon)); // this one fails
		assertThat(exp.log(new BigDecimal(6.091e71, middle)).doubleValue(), isAbsolutely(165.290353874, Epsilon));
		assertThat(exp.log(new BigDecimal(6.091e71, ourDefault)).doubleValue(), isAbsolutely(165.290353874, Epsilon));

		assertNotEquals(exp.log(new BigDecimal(1.245e111, smallest)).doubleValue(), isAbsolutely(255.883339335, Epsilon)); // this one fails
		assertNotEquals(exp.log(new BigDecimal(1.245e111, middle)).doubleValue(), isAbsolutely(255.883339335, Epsilon)); // this one fails
		assertThat(exp.log(new BigDecimal(1.245e111, ourDefault)).doubleValue(), isAbsolutely(255.883339335, 0.1)); // errors are starting to increase!
		assertThat(exp.log(new BigDecimal(1.245e111, new MathContext(75, RoundingMode.HALF_UP))).doubleValue(), isAbsolutely(255.883339335, 0.1)); // errors are starting to increase!

		assertNotEquals(exp.log(new BigDecimal(5.927e300, smallest)).doubleValue(), isAbsolutely(692.555046081, Epsilon)); // this one fails
		assertNotEquals(exp.log(new BigDecimal(5.927e300, middle)).doubleValue(), isAbsolutely(692.555046081, Epsilon)); // this one fails
		assertNotEquals(exp.log(new BigDecimal(5.927e300, ourDefault)).doubleValue(), isAbsolutely(692.555046081, Epsilon)); // this one fails
		assertThat(exp.log(new BigDecimal(5.927e300, new MathContext(264, RoundingMode.HALF_UP))).doubleValue(), isAbsolutely(692.555046081, 0.1));
		assertThat(exp.log(new BigDecimal(5.927e300, new MathContext(300, RoundingMode.HALF_UP))).doubleValue(), isAbsolutely(692.555046081, Epsilon));
	}

	@Test
	public void testLogStable(){
		final double Epsilon = 1e-5;
		final MathContext smallest = new MathContext(8, RoundingMode.HALF_UP);
		final MathContext middle = new MathContext(64, RoundingMode.HALF_UP);
		final MathContext ourDefault = new MathContext(100, RoundingMode.HALF_UP);

		ExpFunction exp = new ExpFunction();
		assertThat(exp.logStable(exp.exp(-10)).doubleValue(), isAbsolutely(-10.0, Epsilon));
		assertThat(exp.logStable(exp.exp(-7.5)).doubleValue(), isAbsolutely(-7.5, Epsilon));
		assertThat(exp.logStable(exp.exp(-5.0)).doubleValue(), isAbsolutely(-5.0, Epsilon));
		assertThat(exp.logStable(exp.exp(-2.5)).doubleValue(), isAbsolutely(-2.5, Epsilon));
		assertThat(exp.logStable(exp.exp(-1.0)).doubleValue(), isAbsolutely(-1.0, Epsilon));
		assertThat(exp.logStable(exp.exp(-0.5)).doubleValue(), isAbsolutely(-0.5, Epsilon));
		assertThat(exp.logStable(exp.exp(0.0)).doubleValue(), isAbsolutely(0.0, Epsilon));
		assertThat(exp.logStable(exp.exp(0.5)).doubleValue(), isAbsolutely(0.5, Epsilon));
		assertThat(exp.logStable(exp.exp(1.0)).doubleValue(), isAbsolutely(1.0, Epsilon));
		assertThat(exp.logStable(exp.exp(2.5)).doubleValue(), isAbsolutely(2.5, Epsilon));
		assertThat(exp.logStable(exp.exp(5.0)).doubleValue(), isAbsolutely(5.0, Epsilon));
		assertThat(exp.logStable(exp.exp(7.5)).doubleValue(), isAbsolutely(7.5, Epsilon));
		assertThat(exp.logStable(exp.exp(10.0)).doubleValue(), isAbsolutely(10.0, Epsilon));

		assertNotEquals(exp.logStable(new BigDecimal(6.091e71, smallest)).doubleValue(), isAbsolutely(165.290353874, Epsilon)); // this one fails
		assertThat(exp.logStable(new BigDecimal(6.091e71, middle)).doubleValue(), isAbsolutely(165.290353874, Epsilon));
		assertThat(exp.logStable(new BigDecimal(6.091e71, ourDefault)).doubleValue(), isAbsolutely(165.290353874, Epsilon));

		assertNotEquals(exp.logStable(new BigDecimal(1.245e111, smallest)).doubleValue(), isAbsolutely(255.883339335, Epsilon)); // this one fails
		assertNotEquals(exp.logStable(new BigDecimal(1.245e111, middle)).doubleValue(), isAbsolutely(255.883339335, Epsilon)); // this one fails
		assertThat(exp.logStable(new BigDecimal(1.245e111, ourDefault)).doubleValue(), isAbsolutely(255.883339335, 0.1)); // errors are starting to increase!
		assertThat(exp.logStable(new BigDecimal(1.245e111, new MathContext(75, RoundingMode.HALF_UP))).doubleValue(), isAbsolutely(255.883339335, 0.1)); // errors are starting to increase!

		assertNotEquals(exp.logStable(new BigDecimal(5.927e300, smallest)).doubleValue(), isAbsolutely(692.555046081, Epsilon)); // this one fails
		assertNotEquals(exp.logStable(new BigDecimal(5.927e300, middle)).doubleValue(), isAbsolutely(692.555046081, Epsilon)); // this one fails
		assertNotEquals(exp.logStable(new BigDecimal(5.927e300, ourDefault)).doubleValue(), isAbsolutely(692.555046081, Epsilon)); // this one fails
		assertThat(exp.logStable(new BigDecimal(5.927e300, new MathContext(264, RoundingMode.HALF_UP))).doubleValue(), isAbsolutely(692.555046081, 0.1));
		assertThat(exp.logStable(new BigDecimal(5.927e300, new MathContext(300, RoundingMode.HALF_UP))).doubleValue(), isAbsolutely(692.555046081, Epsilon));

	}
}
