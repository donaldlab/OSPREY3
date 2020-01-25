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
		final double Epsilon = 1e-8;
		final MathContext smallest = new MathContext(8, RoundingMode.HALF_UP);
		final MathContext ourDefault = new MathContext(64, RoundingMode.HALF_UP);
		final MathContext large = new MathContext(100, RoundingMode.HALF_UP);

		ExpFunction expBad = new ExpFunction();
		ExpFunction expDefault = new ExpFunction();
		ExpFunction expGood = new ExpFunction();

		assertThat(expDefault.log(expDefault.exp(-10)).doubleValue(), isAbsolutely(-10.0, Epsilon));
		assertThat(expDefault.log(expDefault.exp(-7.5)).doubleValue(), isAbsolutely(-7.5, Epsilon));
		assertThat(expDefault.log(expDefault.exp(-5.0)).doubleValue(), isAbsolutely(-5.0, Epsilon));
		assertThat(expDefault.log(expDefault.exp(-2.5)).doubleValue(), isAbsolutely(-2.5, Epsilon));
		assertThat(expDefault.log(expDefault.exp(-1.0)).doubleValue(), isAbsolutely(-1.0, Epsilon));
		assertThat(expDefault.log(expDefault.exp(-0.5)).doubleValue(), isAbsolutely(-0.5, Epsilon));
		assertThat(expDefault.log(expDefault.exp(0.0)).doubleValue(), isAbsolutely(0.0, Epsilon));
		assertThat(expDefault.log(expDefault.exp(0.5)).doubleValue(), isAbsolutely(0.5, Epsilon));
		assertThat(expDefault.log(expDefault.exp(1.0)).doubleValue(), isAbsolutely(1.0, Epsilon));
		assertThat(expDefault.log(expDefault.exp(2.5)).doubleValue(), isAbsolutely(2.5, Epsilon));
		assertThat(expDefault.log(expDefault.exp(5.0)).doubleValue(), isAbsolutely(5.0, Epsilon));
		assertThat(expDefault.log(expDefault.exp(7.5)).doubleValue(), isAbsolutely(7.5, Epsilon));
		assertThat(expDefault.log(expDefault.exp(10.0)).doubleValue(), isAbsolutely(10.0, Epsilon));

		// this log function will fail if the value's exponent is much greater than the precision + 36, empirically
		// testing the effect of precision on the logs
		assertThat(expBad.log(new BigDecimal(6.091e71)).doubleValue(), isAbsolutely(165.290353874, Epsilon));
		assertThat(expDefault.log(new BigDecimal(6.091e71)).doubleValue(), isAbsolutely(165.290353874, Epsilon));
		assertThat(expGood.log(new BigDecimal(6.091e71)).doubleValue(), isAbsolutely(165.290353874, Epsilon));

		assertThat(expBad.log(new BigDecimal(1.245e111)).doubleValue(), isAbsolutely(255.806080852, Epsilon)); // this one fails
		assertThat(expDefault.log(new BigDecimal(1.245e111)).doubleValue(), isAbsolutely(255.806080852, Epsilon)); // this one fails
		assertThat(expGood.log(new BigDecimal(1.245e111)).doubleValue(), isAbsolutely(255.806080852, Epsilon)); // this one fails

		assertThat(expBad.log(new BigDecimal(5.927e300)).doubleValue(), isAbsolutely(692.555046081, Epsilon));
		assertThat(expDefault.log(new BigDecimal(5.927e300)).doubleValue(), isAbsolutely(692.555046081, Epsilon));
		assertThat(expGood.log(new BigDecimal(5.927e300)).doubleValue(), isAbsolutely(692.555046081, Epsilon));
	}

}
