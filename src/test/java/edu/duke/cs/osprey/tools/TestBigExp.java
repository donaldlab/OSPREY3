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

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.*;


import ch.obermuhlner.math.big.BigDecimalMath;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.function.BiFunction;
import java.util.function.Function;


public class TestBigExp {

	private void assertBigExp(BigExp obs, String exp) {
		assertThat(obs.toString(2), is(exp));
	}

	@Test
	public void fromDouble() {

		assertBigExp(new BigExp(Double.POSITIVE_INFINITY), "Infinity");
		assertBigExp(new BigExp(Double.NEGATIVE_INFINITY), "-Infinity");
		assertBigExp(new BigExp(Double.NaN), "NaN");

		assertBigExp(new BigExp(-1e101), "-1.00e101");
		assertBigExp(new BigExp(-1e100), "-1.00e100");
		assertBigExp(new BigExp(-1e99), "-1.00e99");
		assertBigExp(new BigExp(-1e11), "-1.00e11");
		assertBigExp(new BigExp(-1e10), "-1.00e10");
		assertBigExp(new BigExp(-1e9), "-1.00e9");
		assertBigExp(new BigExp(-1e3), "-1.00e3");
		assertBigExp(new BigExp(-100.0), "-1.00e2");
		assertBigExp(new BigExp(-10.0), "-1.00e1");
		assertBigExp(new BigExp(-1.0), "-1.00e0");

		assertBigExp(new BigExp(-0.1), "-1.00e-1");
		assertBigExp(new BigExp(-0.01), "-1.00e-2");
		assertBigExp(new BigExp(-1e-3), "-1.00e-3");
		assertBigExp(new BigExp(-1e-9), "-1.00e-9");
		assertBigExp(new BigExp(-1e-10), "-1.00e-10");
		assertBigExp(new BigExp(-1e-11), "-1.00e-11");
		assertBigExp(new BigExp(-1e-99), "-1.00e-99");
		assertBigExp(new BigExp(-1e-100), "-1.00e-100");
		assertBigExp(new BigExp(-1e-101), "-1.00e-101");

		assertBigExp(new BigExp(0.0), "0.00e0");

		assertBigExp(new BigExp(1e-101), "1.00e-101");
		assertBigExp(new BigExp(1e-100), "1.00e-100");
		assertBigExp(new BigExp(1e-99), "1.00e-99");
		assertBigExp(new BigExp(1e-11), "1.00e-11");
		assertBigExp(new BigExp(1e-10), "1.00e-10");
		assertBigExp(new BigExp(1e-9), "1.00e-9");
		assertBigExp(new BigExp(1e-3), "1.00e-3");
		assertBigExp(new BigExp(0.01), "1.00e-2");
		assertBigExp(new BigExp(0.1), "1.00e-1");

		assertBigExp(new BigExp(1.0), "1.00e0");
		assertBigExp(new BigExp(10.0), "1.00e1");
		assertBigExp(new BigExp(100.0), "1.00e2");
		assertBigExp(new BigExp(1e3), "1.00e3");
		assertBigExp(new BigExp(1e9), "1.00e9");
		assertBigExp(new BigExp(1e10), "1.00e10");
		assertBigExp(new BigExp(1e11), "1.00e11");
		assertBigExp(new BigExp(1e99), "1.00e99");
		assertBigExp(new BigExp(1e100), "1.00e100");
		assertBigExp(new BigExp(1e101), "1.00e101");

		assertBigExp(new BigExp(4.23), "4.23e0");
		assertBigExp(new BigExp(9.94e302), "9.94e302");
		assertBigExp(new BigExp(9.82e-301), "9.82e-301");
		assertBigExp(new BigExp(-9.94e302), "-9.94e302");
		assertBigExp(new BigExp(-9.82e-301), "-9.82e-301");
	}

	@Test
	public void fromBigDecimal() {

		assertBigExp(new BigExp(MathTools.BigPositiveInfinity), "Infinity");
		assertBigExp(new BigExp(MathTools.BigNegativeInfinity), "-Infinity");
		assertBigExp(new BigExp(MathTools.BigNaN), "NaN");

		assertBigExp(new BigExp(BigDecimal.ZERO), "0.00e0");
		assertBigExp(new BigExp(BigDecimal.ONE), "1.00e0");
		assertBigExp(new BigExp(BigDecimal.TEN), "1.00e1");

		assertBigExp(new BigExp(new BigDecimal("1.23e10")), "1.23e10");
		assertBigExp(new BigExp(new BigDecimal("1.23e-10")), "1.23e-10");
		assertBigExp(new BigExp(new BigDecimal("-1.23e10")), "-1.23e10");
		assertBigExp(new BigExp(new BigDecimal("-1.23e-10")), "-1.23e-10");

		assertBigExp(new BigExp(new BigDecimal("9.58e400")), "9.58e400");
		assertBigExp(new BigExp(new BigDecimal("9.58e-400")), "9.58e-400");
		assertBigExp(new BigExp(new BigDecimal("-9.58e400")), "-9.58e400");
		assertBigExp(new BigExp(new BigDecimal("-9.58e-400")), "-9.58e-400");
	}

	@Test
	public void toBigDecimal() {

		assertThat(new BigExp(Double.POSITIVE_INFINITY).toBigDecimal(), is(MathTools.BigPositiveInfinity));
		assertThat(new BigExp(Double.NEGATIVE_INFINITY).toBigDecimal(), is(MathTools.BigNegativeInfinity));
		assertThat(new BigExp(Double.NaN).toBigDecimal(), is(MathTools.BigNaN));

		Function<BigExp,String> f = bigexp -> String.format("%.2e", bigexp.toBigDecimal());

		assertThat(f.apply(new BigExp(1.23, -1000)), is("1.23e-1000"));
		assertThat(f.apply(new BigExp(1.23, -100)), is("1.23e-100"));
		assertThat(f.apply(new BigExp(1.23, -10)), is("1.23e-10"));
		assertThat(f.apply(new BigExp(1.23, -1)), is("1.23e-01"));
		assertThat(f.apply(new BigExp(1.23, 0)), is("1.23e+00"));
		assertThat(f.apply(new BigExp(1.23, 1)), is("1.23e+01"));
		assertThat(f.apply(new BigExp(1.23, 10)), is("1.23e+10"));
		assertThat(f.apply(new BigExp(1.23, 100)), is("1.23e+100"));
		assertThat(f.apply(new BigExp(1.23, 1000)), is("1.23e+1000"));

		assertThat(f.apply(new BigExp(-1.23, -1000)), is("-1.23e-1000"));
		assertThat(f.apply(new BigExp(-1.23, -100)), is("-1.23e-100"));
		assertThat(f.apply(new BigExp(-1.23, -10)), is("-1.23e-10"));
		assertThat(f.apply(new BigExp(-1.23, -1)), is("-1.23e-01"));
		assertThat(f.apply(new BigExp(-1.23, 0)), is("-1.23e+00"));
		assertThat(f.apply(new BigExp(-1.23, 1)), is("-1.23e+01"));
		assertThat(f.apply(new BigExp(-1.23, 10)), is("-1.23e+10"));
		assertThat(f.apply(new BigExp(-1.23, 100)), is("-1.23e+100"));
		assertThat(f.apply(new BigExp(-1.23, 1000)), is("-1.23e+1000"));
	}

	@Test
	public void fromBigDecimalPrecise() {

		// check accuracy to 16 significant digits (that's about all double can represent)
		Function<Double,String> d = val -> String.format("%.15e", val);
		Function<BigExp,String> be = val -> String.format("%.15e", val.toBigDecimal());

		assertThat(d.apply(1.1), is("1.100000000000000e+00"));
		assertThat(be.apply(new BigExp(new BigDecimal("1.1e0"))), is("1.100000000000000e+00"));

		assertThat(d.apply(1.111111), is("1.111111000000000e+00"));
		assertThat(be.apply(new BigExp(new BigDecimal("1.111111e0"))), is("1.111111000000000e+00"));

		assertThat(d.apply(1.11111111111), is("1.111111111110000e+00"));
		assertThat(be.apply(new BigExp(new BigDecimal("1.11111111111e0"))), is("1.111111111110000e+00"));

		assertThat(d.apply(1.111111111111111), is("1.111111111111111e+00"));
		assertThat(be.apply(new BigExp(new BigDecimal("1.111111111111111e0"))), is("1.111111111111111e+00"));

		assertThat(d.apply(1.111111111111111e1), is("1.111111111111111e+01"));
		assertThat(be.apply(new BigExp(new BigDecimal("1.111111111111111e1"))), is("1.111111111111111e+01"));

		assertThat(d.apply(1.111111111111111e100), is("1.111111111111111e+100"));
		assertThat(be.apply(new BigExp(new BigDecimal("1.111111111111111e100"))), is("1.111111111111111e+100"));

		assertThat(d.apply(1.111111111111111e299), is("1.111111111111111e+299"));
		assertThat(be.apply(new BigExp(new BigDecimal("1.111111111111111e299"))), is("1.111111111111111e+299"));

		assertThat(d.apply(1.111111111111111e-299), is("1.111111111111111e-299"));
		assertThat(be.apply(new BigExp(new BigDecimal("1.111111111111111e-299"))), is("1.111111111111111e-299"));
	}

	@Test
	public void toBigDecimalPrecise() {

		// check accuracy to 16 significant digits (that's about all double can represent)
		Function<Double,String> d = val -> String.format("%.15e", val);
		Function<BigExp,String> bd = val -> String.format("%.15e", val.toBigDecimal());

		assertThat(d.apply(1.1), is("1.100000000000000e+00"));
		assertThat(bd.apply(new BigExp(1.1, 0)), is("1.100000000000000e+00"));

		assertThat(d.apply(1.111111), is("1.111111000000000e+00"));
		assertThat(bd.apply(new BigExp(1.111111, 0)), is("1.111111000000000e+00"));

		assertThat(d.apply(1.11111111111), is("1.111111111110000e+00"));
		assertThat(bd.apply(new BigExp(1.11111111111, 0)), is("1.111111111110000e+00"));

		assertThat(d.apply(1.111111111111111), is("1.111111111111111e+00"));
		assertThat(bd.apply(new BigExp(1.111111111111111, 0)), is("1.111111111111111e+00"));

		assertThat(d.apply(1.111111111111111e1), is("1.111111111111111e+01"));
		assertThat(bd.apply(new BigExp(1.111111111111111, 1)), is("1.111111111111111e+01"));

		assertThat(d.apply(1.111111111111111e100), is("1.111111111111111e+100"));
		assertThat(bd.apply(new BigExp(1.111111111111111, 100)), is("1.111111111111111e+100"));

		assertThat(d.apply(1.111111111111111e299), is("1.111111111111111e+299"));
		assertThat(bd.apply(new BigExp(1.111111111111111, 299)), is("1.111111111111111e+299"));

		assertThat(d.apply(1.111111111111111e-299), is("1.111111111111111e-299"));
		assertThat(bd.apply(new BigExp(1.111111111111111, -299)), is("1.111111111111111e-299"));
	}

	@Test
	public void toBigDecimalUpconvert() {

		// check accuracy to 20 significant digits (that's a bit more than double can represent)
		// but make sure the extra precision is all zeros
		MathContext mathContext = new MathContext(20, RoundingMode.HALF_UP);
		Function<Double,String> d = val -> String.format("%.19e", val);
		Function<BigExp,String> bd = val -> String.format("%.19e", val.toBigDecimal(mathContext));

		assertThat(d.apply(1.1), is("1.1000000000000000000e+00"));
		assertThat(bd.apply(new BigExp(1.1, 0)), is("1.1000000000000000000e+00"));

		assertThat(d.apply(1.111111), is("1.1111110000000000000e+00"));
		assertThat(bd.apply(new BigExp(1.111111, 0)), is("1.1111110000000000000e+00"));

		assertThat(d.apply(1.11111111111), is("1.1111111111100000000e+00"));
		assertThat(bd.apply(new BigExp(1.11111111111, 0)), is("1.1111111111100000000e+00"));

		assertThat(d.apply(1.111111111111111), is("1.1111111111111110000e+00"));
		assertThat(bd.apply(new BigExp(1.111111111111111, 0)), is("1.1111111111111110000e+00"));

		assertThat(d.apply(1.111111111111111e1), is("1.1111111111111110000e+01"));
		assertThat(bd.apply(new BigExp(1.111111111111111, 1)), is("1.1111111111111110000e+01"));

		assertThat(d.apply(1.111111111111111e100), is("1.1111111111111110000e+100"));
		assertThat(bd.apply(new BigExp(1.111111111111111, 100)), is("1.1111111111111110000e+100"));

		assertThat(d.apply(1.111111111111111e299), is("1.1111111111111110000e+299"));
		assertThat(bd.apply(new BigExp(1.111111111111111, 299)), is("1.1111111111111110000e+299"));

		assertThat(d.apply(1.111111111111111e-299), is("1.1111111111111110000e-299"));
		assertThat(bd.apply(new BigExp(1.111111111111111, -299)), is("1.1111111111111110000e-299"));
	}

	@Test
	public void mult() {

		BiFunction<BigExp,BigExp,BigExp> mult = (a, b) -> {
			BigExp out = new BigExp(a);
			out.mult(b);
			return out;
		};

		assertBigExp(mult.apply(new BigExp(0.0, 0), new BigExp(0.0, 0)), "0.00e0");

		assertBigExp(mult.apply(new BigExp(1.0, 0), new BigExp(1.0, 0)), "1.00e0");

		assertBigExp(mult.apply(new BigExp(3.29, -1000), new BigExp(0.0, 0)), "0.00e0");
		assertBigExp(mult.apply(new BigExp(3.29, -100), new BigExp(0.0, 0)), "0.00e0");
		assertBigExp(mult.apply(new BigExp(3.29, -10), new BigExp(0.0, 0)), "0.00e0");
		assertBigExp(mult.apply(new BigExp(3.29, -1), new BigExp(0.0, 0)), "0.00e0");
		assertBigExp(mult.apply(new BigExp(3.29, 0), new BigExp(0.0, 0)), "0.00e0");
		assertBigExp(mult.apply(new BigExp(3.29, 1), new BigExp(0.0, 0)), "0.00e0");
		assertBigExp(mult.apply(new BigExp(3.29, 10), new BigExp(0.0, 0)), "0.00e0");
		assertBigExp(mult.apply(new BigExp(3.29, 100), new BigExp(0.0, 0)), "0.00e0");
		assertBigExp(mult.apply(new BigExp(3.29, 1000), new BigExp(0.0, 0)), "0.00e0");

		assertBigExp(mult.apply(new BigExp(-3.29, -1000), new BigExp(0.0, 0)), "-0.00e0");
		assertBigExp(mult.apply(new BigExp(-3.29, -100), new BigExp(0.0, 0)), "-0.00e0");
		assertBigExp(mult.apply(new BigExp(-3.29, -10), new BigExp(0.0, 0)), "-0.00e0");
		assertBigExp(mult.apply(new BigExp(-3.29, -1), new BigExp(0.0, 0)), "-0.00e0");
		assertBigExp(mult.apply(new BigExp(-3.29, 0), new BigExp(0.0, 0)), "-0.00e0");
		assertBigExp(mult.apply(new BigExp(-3.29, 1), new BigExp(0.0, 0)), "-0.00e0");
		assertBigExp(mult.apply(new BigExp(-3.29, 10), new BigExp(0.0, 0)), "-0.00e0");
		assertBigExp(mult.apply(new BigExp(-3.29, 100), new BigExp(0.0, 0)), "-0.00e0");
		assertBigExp(mult.apply(new BigExp(-3.29, 1000), new BigExp(0.0, 0)), "-0.00e0");


		assertBigExp(mult.apply(new BigExp(8.43, -1000), new BigExp(3.92, -1000)), "3.30e-1999");
		assertBigExp(mult.apply(new BigExp(8.43, -100), new BigExp(3.92, -1000)), "3.30e-1099");
		assertBigExp(mult.apply(new BigExp(8.43, -10), new BigExp(3.92, -1000)), "3.30e-1009");
		assertBigExp(mult.apply(new BigExp(8.43, -1), new BigExp(3.92, -1000)), "3.30e-1000");
		assertBigExp(mult.apply(new BigExp(8.43, 0), new BigExp(3.92, -1000)), "3.30e-999");
		assertBigExp(mult.apply(new BigExp(8.43, 1), new BigExp(3.92, -1000)), "3.30e-998");
		assertBigExp(mult.apply(new BigExp(8.43, 10), new BigExp(3.92, -1000)), "3.30e-989");
		assertBigExp(mult.apply(new BigExp(8.43, 100), new BigExp(3.92, -1000)), "3.30e-899");
		assertBigExp(mult.apply(new BigExp(8.43, 1000), new BigExp(3.92, -1000)), "3.30e1");

		assertBigExp(mult.apply(new BigExp(8.43, -1000), new BigExp(3.92, -100)), "3.30e-1099");
		assertBigExp(mult.apply(new BigExp(8.43, -100), new BigExp(3.92, -100)), "3.30e-199");
		assertBigExp(mult.apply(new BigExp(8.43, -10), new BigExp(3.92, -100)), "3.30e-109");
		assertBigExp(mult.apply(new BigExp(8.43, -1), new BigExp(3.92, -100)), "3.30e-100");
		assertBigExp(mult.apply(new BigExp(8.43, 0), new BigExp(3.92, -100)), "3.30e-99");
		assertBigExp(mult.apply(new BigExp(8.43, 1), new BigExp(3.92, -100)), "3.30e-98");
		assertBigExp(mult.apply(new BigExp(8.43, 10), new BigExp(3.92, -100)), "3.30e-89");
		assertBigExp(mult.apply(new BigExp(8.43, 100), new BigExp(3.92, -100)), "3.30e1");
		assertBigExp(mult.apply(new BigExp(8.43, 1000), new BigExp(3.92, -100)), "3.30e901");

		assertBigExp(mult.apply(new BigExp(8.43, -1000), new BigExp(3.92, -10)), "3.30e-1009");
		assertBigExp(mult.apply(new BigExp(8.43, -100), new BigExp(3.92, -10)), "3.30e-109");
		assertBigExp(mult.apply(new BigExp(8.43, -10), new BigExp(3.92, -10)), "3.30e-19");
		assertBigExp(mult.apply(new BigExp(8.43, -1), new BigExp(3.92, -10)), "3.30e-10");
		assertBigExp(mult.apply(new BigExp(8.43, 0), new BigExp(3.92, -10)), "3.30e-9");
		assertBigExp(mult.apply(new BigExp(8.43, 1), new BigExp(3.92, -10)), "3.30e-8");
		assertBigExp(mult.apply(new BigExp(8.43, 10), new BigExp(3.92, -10)), "3.30e1");
		assertBigExp(mult.apply(new BigExp(8.43, 100), new BigExp(3.92, -10)), "3.30e91");
		assertBigExp(mult.apply(new BigExp(8.43, 1000), new BigExp(3.92, -10)), "3.30e991");

		assertBigExp(mult.apply(new BigExp(8.43, -1000), new BigExp(3.92, -1)), "3.30e-1000");
		assertBigExp(mult.apply(new BigExp(8.43, -100), new BigExp(3.92, -1)), "3.30e-100");
		assertBigExp(mult.apply(new BigExp(8.43, -10), new BigExp(3.92, -1)), "3.30e-10");
		assertBigExp(mult.apply(new BigExp(8.43, -1), new BigExp(3.92, -1)), "3.30e-1");
		assertBigExp(mult.apply(new BigExp(8.43, 0), new BigExp(3.92, -1)), "3.30e0");
		assertBigExp(mult.apply(new BigExp(8.43, 1), new BigExp(3.92, -1)), "3.30e1");
		assertBigExp(mult.apply(new BigExp(8.43, 10), new BigExp(3.92, -1)), "3.30e10");
		assertBigExp(mult.apply(new BigExp(8.43, 100), new BigExp(3.92, -1)), "3.30e100");
		assertBigExp(mult.apply(new BigExp(8.43, 1000), new BigExp(3.92, -1)), "3.30e1000");

		assertBigExp(mult.apply(new BigExp(8.43, -1000), new BigExp(3.92, 0)), "3.30e-999");
		assertBigExp(mult.apply(new BigExp(8.43, -100), new BigExp(3.92, 0)), "3.30e-99");
		assertBigExp(mult.apply(new BigExp(8.43, -10), new BigExp(3.92, 0)), "3.30e-9");
		assertBigExp(mult.apply(new BigExp(8.43, -1), new BigExp(3.92, 0)), "3.30e0");
		assertBigExp(mult.apply(new BigExp(8.43, 0), new BigExp(3.92, 0)), "3.30e1");
		assertBigExp(mult.apply(new BigExp(8.43, 1), new BigExp(3.92, 0)), "3.30e2");
		assertBigExp(mult.apply(new BigExp(8.43, 10), new BigExp(3.92, 0)), "3.30e11");
		assertBigExp(mult.apply(new BigExp(8.43, 100), new BigExp(3.92, 0)), "3.30e101");
		assertBigExp(mult.apply(new BigExp(8.43, 1000), new BigExp(3.92, 0)), "3.30e1001");

		assertBigExp(mult.apply(new BigExp(8.43, -1000), new BigExp(3.92, 1)), "3.30e-998");
		assertBigExp(mult.apply(new BigExp(8.43, -100), new BigExp(3.92, 1)), "3.30e-98");
		assertBigExp(mult.apply(new BigExp(8.43, -10), new BigExp(3.92, 1)), "3.30e-8");
		assertBigExp(mult.apply(new BigExp(8.43, -1), new BigExp(3.92, 1)), "3.30e1");
		assertBigExp(mult.apply(new BigExp(8.43, 0), new BigExp(3.92, 1)), "3.30e2");
		assertBigExp(mult.apply(new BigExp(8.43, 1), new BigExp(3.92, 1)), "3.30e3");
		assertBigExp(mult.apply(new BigExp(8.43, 10), new BigExp(3.92, 1)), "3.30e12");
		assertBigExp(mult.apply(new BigExp(8.43, 100), new BigExp(3.92, 1)), "3.30e102");
		assertBigExp(mult.apply(new BigExp(8.43, 1000), new BigExp(3.92, 1)), "3.30e1002");

		assertBigExp(mult.apply(new BigExp(8.43, -1000), new BigExp(3.92, 10)), "3.30e-989");
		assertBigExp(mult.apply(new BigExp(8.43, -100), new BigExp(3.92, 10)), "3.30e-89");
		assertBigExp(mult.apply(new BigExp(8.43, -10), new BigExp(3.92, 10)), "3.30e1");
		assertBigExp(mult.apply(new BigExp(8.43, -1), new BigExp(3.92, 10)), "3.30e10");
		assertBigExp(mult.apply(new BigExp(8.43, 0), new BigExp(3.92, 10)), "3.30e11");
		assertBigExp(mult.apply(new BigExp(8.43, 1), new BigExp(3.92, 10)), "3.30e12");
		assertBigExp(mult.apply(new BigExp(8.43, 10), new BigExp(3.92, 10)), "3.30e21");
		assertBigExp(mult.apply(new BigExp(8.43, 100), new BigExp(3.92, 10)), "3.30e111");
		assertBigExp(mult.apply(new BigExp(8.43, 1000), new BigExp(3.92, 10)), "3.30e1011");

		assertBigExp(mult.apply(new BigExp(8.43, -1000), new BigExp(3.92, 100)), "3.30e-899");
		assertBigExp(mult.apply(new BigExp(8.43, -100), new BigExp(3.92, 100)), "3.30e1");
		assertBigExp(mult.apply(new BigExp(8.43, -10), new BigExp(3.92, 100)), "3.30e91");
		assertBigExp(mult.apply(new BigExp(8.43, -1), new BigExp(3.92, 100)), "3.30e100");
		assertBigExp(mult.apply(new BigExp(8.43, 0), new BigExp(3.92, 100)), "3.30e101");
		assertBigExp(mult.apply(new BigExp(8.43, 1), new BigExp(3.92, 100)), "3.30e102");
		assertBigExp(mult.apply(new BigExp(8.43, 10), new BigExp(3.92, 100)), "3.30e111");
		assertBigExp(mult.apply(new BigExp(8.43, 100), new BigExp(3.92, 100)), "3.30e201");
		assertBigExp(mult.apply(new BigExp(8.43, 1000), new BigExp(3.92, 100)), "3.30e1101");

		assertBigExp(mult.apply(new BigExp(8.43, -1000), new BigExp(3.92, 1000)), "3.30e1");
		assertBigExp(mult.apply(new BigExp(8.43, -100), new BigExp(3.92, 1000)), "3.30e901");
		assertBigExp(mult.apply(new BigExp(8.43, -10), new BigExp(3.92, 1000)), "3.30e991");
		assertBigExp(mult.apply(new BigExp(8.43, -1), new BigExp(3.92, 1000)), "3.30e1000");
		assertBigExp(mult.apply(new BigExp(8.43, 0), new BigExp(3.92, 1000)), "3.30e1001");
		assertBigExp(mult.apply(new BigExp(8.43, 1), new BigExp(3.92, 1000)), "3.30e1002");
		assertBigExp(mult.apply(new BigExp(8.43, 10), new BigExp(3.92, 1000)), "3.30e1011");
		assertBigExp(mult.apply(new BigExp(8.43, 100), new BigExp(3.92, 1000)), "3.30e1101");
		assertBigExp(mult.apply(new BigExp(8.43, 1000), new BigExp(3.92, 1000)), "3.30e2001");

		assertThat(8.43e200 * 3.92e200, is(Double.POSITIVE_INFINITY));
		assertBigExp(mult.apply(new BigExp(8.43, 200), new BigExp(3.92, 200)), "3.30e401");
		assertBigExp(mult.apply(new BigExp(8.43e200), new BigExp(3.92e200)), "3.30e401");

		assertThat(-8.43e200 * 3.92e200, is(Double.NEGATIVE_INFINITY));
		assertBigExp(mult.apply(new BigExp(-8.43, 200), new BigExp(3.92, 200)), "-3.30e401");
		assertBigExp(mult.apply(new BigExp(-8.43e200), new BigExp(3.92e200)), "-3.30e401");

		assertThat(8.43e-200 * 3.92e-200, is(0.0));
		assertBigExp(mult.apply(new BigExp(8.43, -200), new BigExp(3.92, -200)), "3.30e-399");
		assertBigExp(mult.apply(new BigExp(8.43e-200), new BigExp(3.92e-200)), "3.30e-399");

		assertThat(-8.43e-200 * 3.92e-200, is(-0.0));
		assertBigExp(mult.apply(new BigExp(-8.43, -200), new BigExp(3.92, -200)), "-3.30e-399");
		assertBigExp(mult.apply(new BigExp(-8.43e-200), new BigExp(3.92e-200)), "-3.30e-399");
	}

	@Test
	public void compare() {

		assertThat(new BigExp(0.0, 0).compareTo(new BigExp(0.0, 0)), is(0));

		assertThat(new BigExp(6.66, -100).compareTo(new BigExp(0.0, 0)), is(+1));
		assertThat(new BigExp(6.66, -10).compareTo(new BigExp(0.0, 0)), is(+1));
		assertThat(new BigExp(6.66, -1).compareTo(new BigExp(0.0, 0)), is(+1));
		assertThat(new BigExp(6.66,  0).compareTo(new BigExp(0.0, 0)), is(+1));
		assertThat(new BigExp(6.66, +1).compareTo(new BigExp(0.0, 0)), is(+1));
		assertThat(new BigExp(6.66, +10).compareTo(new BigExp(0.0, 0)), is(+1));
		assertThat(new BigExp(6.66, +100).compareTo(new BigExp(0.0, 0)), is(+1));

		assertThat(new BigExp(-6.66, -100).compareTo(new BigExp(0.0, 0)), is(-1));
		assertThat(new BigExp(-6.66, -10).compareTo(new BigExp(0.0, 0)), is(-1));
		assertThat(new BigExp(-6.66, -1).compareTo(new BigExp(0.0, 0)), is(-1));
		assertThat(new BigExp(-6.66,  0).compareTo(new BigExp(0.0, 0)), is(-1));
		assertThat(new BigExp(-6.66, +1).compareTo(new BigExp(0.0, 0)), is(-1));
		assertThat(new BigExp(-6.66, +10).compareTo(new BigExp(0.0, 0)), is(-1));
		assertThat(new BigExp(-6.66, +100).compareTo(new BigExp(0.0, 0)), is(-1));

		assertThat(new BigExp(4.93, 0).compareTo(new BigExp(2.64, 0)), is(+1));
		assertThat(new BigExp(2.64, 0).compareTo(new BigExp(4.93, 0)), is(-1));

		assertThat(new BigExp(4.93e10, 0).compareTo(new BigExp(2.64e10, 0)), is(+1));
		assertThat(new BigExp(2.64e10, 0).compareTo(new BigExp(4.93e10, 0)), is(-1));

		assertThat(new BigExp(4.93e10, 10).compareTo(new BigExp(2.64e10, 5)), is(+1));
		assertThat(new BigExp(2.64e10, 5).compareTo(new BigExp(4.93e10, 10)), is(-1));

		assertThat(new BigExp(4.93e10, 5).compareTo(new BigExp(2.64e10, 10)), is(-1));
		assertThat(new BigExp(2.64e10, 10).compareTo(new BigExp(4.93e10, 5)), is(+1));

		assertThat(new BigExp(4.93e20, 5).compareTo(new BigExp(2.64e10, 10)), is(+1));
		assertThat(new BigExp(2.64e10, 10).compareTo(new BigExp(4.93e20, 5)), is(-1));

		assertThat(new BigExp(4.93e5, 5).compareTo(new BigExp(2.64e5, 6)), is(-1));
		assertThat(new BigExp(2.64e5, 6).compareTo(new BigExp(4.93e5, 5)), is(+1));

		assertThat(new BigExp(4.93, 0).compareTo(new BigExp(2.64, 1)), is(-1));
		assertThat(new BigExp(4.93, 1).compareTo(new BigExp(2.64, 0)), is(+1));
		assertThat(new BigExp(4.93, 1).compareTo(new BigExp(2.64, 1)), is(+1));

		assertThat(new BigExp(4.93,  0).compareTo(new BigExp(2.64, -1)), is(+1));
		assertThat(new BigExp(4.93, -1).compareTo(new BigExp(2.64,  0)), is(-1));
		assertThat(new BigExp(4.93, -1).compareTo(new BigExp(2.64, -1)), is(+1));

		assertThat(new BigExp(6.66, -1).compareTo(new BigExp(6.66, -1)), is( 0));
		assertThat(new BigExp(6.66, -1).compareTo(new BigExp(6.66,  0)), is(-1));
		assertThat(new BigExp(6.66, -1).compareTo(new BigExp(6.66, +1)), is(-1));
		assertThat(new BigExp(6.66,  0).compareTo(new BigExp(6.66, -1)), is(+1));
		assertThat(new BigExp(6.66,  0).compareTo(new BigExp(6.66,  0)), is( 0));
		assertThat(new BigExp(6.66,  0).compareTo(new BigExp(6.66, +1)), is(-1));
		assertThat(new BigExp(6.66, +1).compareTo(new BigExp(6.66, -1)), is(+1));
		assertThat(new BigExp(6.66, +1).compareTo(new BigExp(6.66,  0)), is(+1));
		assertThat(new BigExp(6.66, +1).compareTo(new BigExp(6.66, +1)), is( 0));

		// errors seen in the wild
		assertThat(new BigExp(1.673007, -78675272).compareTo(new BigExp(0.0, 0)), is(+1));
		assertThat(new BigExp(1.851726, 47).compareTo(new BigExp(7.582769, 48)), is(-1));

		// infinity errors
		assertThat(new BigExp(10,1000000).compareTo(new BigExp(0.0, 0)), is(+1));
		assertThat(new BigExp(Double.POSITIVE_INFINITY).compareTo(new BigExp(0.0, 0)), is(+1));
		assertThat(new BigExp(5.852204485365044,7).compareTo(new BigExp(Double.POSITIVE_INFINITY)), is(-1));
		assertThat(new BigExp(0.0,0).compareTo(new BigExp(10,1000000)), is(-1));

	}

	@Test
	public void greaterThan() {

		final double epsilon = 1e-1;

		assertThat(new BigExp(0.0, 0).greaterThan(new BigExp(0.0, 0), epsilon), is(false));

		assertThat(new BigExp(6.66, -100).greaterThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(6.66, -10).greaterThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(6.66, -1).greaterThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(6.66,  0).greaterThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(6.66, +1).greaterThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(6.66, +10).greaterThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(6.66, +100).greaterThan(new BigExp(0.0, 0), epsilon), is(true));

		assertThat(new BigExp(-6.66, -100).greaterThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(-6.66, -10).greaterThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(-6.66, -1).greaterThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(-6.66,  0).greaterThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(-6.66, +1).greaterThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(-6.66, +10).greaterThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(-6.66, +100).greaterThan(new BigExp(0.0, 0), epsilon), is(false));

		assertThat(new BigExp(4.93, 0).greaterThan(new BigExp(2.64, 0), epsilon), is(true));
		assertThat(new BigExp(2.64, 0).greaterThan(new BigExp(4.93, 0), epsilon), is(false));

		assertThat(new BigExp(4.93e10, 0).greaterThan(new BigExp(2.64e10, 0), epsilon), is(true));
		assertThat(new BigExp(2.64e10, 0).greaterThan(new BigExp(4.93e10, 0), epsilon), is(false));

		assertThat(new BigExp(4.93e10, 10).greaterThan(new BigExp(2.64e10, 5), epsilon), is(true));
		assertThat(new BigExp(2.64e10, 5).greaterThan(new BigExp(4.93e10, 10), epsilon), is(false));

		assertThat(new BigExp(4.93e10, 5).greaterThan(new BigExp(2.64e10, 10), epsilon), is(false));
		assertThat(new BigExp(2.64e10, 10).greaterThan(new BigExp(4.93e10, 5), epsilon), is(true));

		assertThat(new BigExp(4.93e20, 5).greaterThan(new BigExp(2.64e10, 10), epsilon), is(true));
		assertThat(new BigExp(2.64e10, 10).greaterThan(new BigExp(4.93e20, 5), epsilon), is(false));

		assertThat(new BigExp(4.93e5, 5).greaterThan(new BigExp(2.64e5, 6), epsilon), is(false));
		assertThat(new BigExp(2.64e5, 6).greaterThan(new BigExp(4.93e5, 5), epsilon), is(true));

		assertThat(new BigExp(4.93, 0).greaterThan(new BigExp(2.64, 1), epsilon), is(false));
		assertThat(new BigExp(4.93, 1).greaterThan(new BigExp(2.64, 0), epsilon), is(true));
		assertThat(new BigExp(4.93, 1).greaterThan(new BigExp(2.64, 1), epsilon), is(true));

		assertThat(new BigExp(4.93,  0).greaterThan(new BigExp(2.64, -1), epsilon), is(true));
		assertThat(new BigExp(4.93, -1).greaterThan(new BigExp(2.64,  0), epsilon), is(false));
		assertThat(new BigExp(4.93, -1).greaterThan(new BigExp(2.64, -1), epsilon), is(true));

		assertThat(new BigExp(6.66, -1).greaterThan(new BigExp(6.66, -1), epsilon), is(false));
		assertThat(new BigExp(6.66, -1).greaterThan(new BigExp(6.66,  0), epsilon), is(false));
		assertThat(new BigExp(6.66, -1).greaterThan(new BigExp(6.66, +1), epsilon), is(false));
		assertThat(new BigExp(6.66,  0).greaterThan(new BigExp(6.66, -1), epsilon), is(true));
		assertThat(new BigExp(6.66,  0).greaterThan(new BigExp(6.66,  0), epsilon), is(false));
		assertThat(new BigExp(6.66,  0).greaterThan(new BigExp(6.66, +1), epsilon), is(false));
		assertThat(new BigExp(6.66, +1).greaterThan(new BigExp(6.66, -1), epsilon), is(true));
		assertThat(new BigExp(6.66, +1).greaterThan(new BigExp(6.66,  0), epsilon), is(true));
		assertThat(new BigExp(6.66, +1).greaterThan(new BigExp(6.66, +1), epsilon), is(false));

		// issues caused by setting the exponent of infinity to zero when normalizing
		assertThat(new BigExp(6.66, +1).greaterThan(new BigExp(Double.NEGATIVE_INFINITY)), is(true));
		assertThat(new BigExp(-6.66, +1).greaterThan(new BigExp(Double.NEGATIVE_INFINITY)), is(true));
	}

	@Test
	public void lessThan() {

		final double epsilon = 1e-1;

		assertThat(new BigExp(0.0, 0).lessThan(new BigExp(0.0, 0), epsilon), is(false));

		assertThat(new BigExp(6.66, -100).lessThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(6.66, -10).lessThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(6.66, -1).lessThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(6.66,  0).lessThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(6.66, +1).lessThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(6.66, +10).lessThan(new BigExp(0.0, 0), epsilon), is(false));
		assertThat(new BigExp(6.66, +100).lessThan(new BigExp(0.0, 0), epsilon), is(false));

		assertThat(new BigExp(-6.66, -100).lessThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(-6.66, -10).lessThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(-6.66, -1).lessThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(-6.66,  0).lessThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(-6.66, +1).lessThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(-6.66, +10).lessThan(new BigExp(0.0, 0), epsilon), is(true));
		assertThat(new BigExp(-6.66, +100).lessThan(new BigExp(0.0, 0), epsilon), is(true));

		assertThat(new BigExp(4.93, 0).lessThan(new BigExp(2.64, 0), epsilon), is(false));
		assertThat(new BigExp(2.64, 0).lessThan(new BigExp(4.93, 0), epsilon), is(true));

		assertThat(new BigExp(4.93e10, 0).lessThan(new BigExp(2.64e10, 0), epsilon), is(false));
		assertThat(new BigExp(2.64e10, 0).lessThan(new BigExp(4.93e10, 0), epsilon), is(true));

		assertThat(new BigExp(4.93e10, 10).lessThan(new BigExp(2.64e10, 5), epsilon), is(false));
		assertThat(new BigExp(2.64e10, 5).lessThan(new BigExp(4.93e10, 10), epsilon), is(true));

		assertThat(new BigExp(4.93e10, 5).lessThan(new BigExp(2.64e10, 10), epsilon), is(true));
		assertThat(new BigExp(2.64e10, 10).lessThan(new BigExp(4.93e10, 5), epsilon), is(false));

		assertThat(new BigExp(4.93e20, 5).lessThan(new BigExp(2.64e10, 10), epsilon), is(false));
		assertThat(new BigExp(2.64e10, 10).lessThan(new BigExp(4.93e20, 5), epsilon), is(true));

		assertThat(new BigExp(4.93e5, 5).lessThan(new BigExp(2.64e5, 6), epsilon), is(true));
		assertThat(new BigExp(2.64e5, 6).lessThan(new BigExp(4.93e5, 5), epsilon), is(false));

		assertThat(new BigExp(4.93, 0).lessThan(new BigExp(2.64, 1), epsilon), is(true));
		assertThat(new BigExp(4.93, 1).lessThan(new BigExp(2.64, 0), epsilon), is(false));
		assertThat(new BigExp(4.93, 1).lessThan(new BigExp(2.64, 1), epsilon), is(false));

		assertThat(new BigExp(4.93,  0).lessThan(new BigExp(2.64, -1), epsilon), is(false));
		assertThat(new BigExp(4.93, -1).lessThan(new BigExp(2.64,  0), epsilon), is(true));
		assertThat(new BigExp(4.93, -1).lessThan(new BigExp(2.64, -1), epsilon), is(false));

		assertThat(new BigExp(6.66, -1).lessThan(new BigExp(6.66, -1), epsilon), is(false));
		assertThat(new BigExp(6.66, -1).lessThan(new BigExp(6.66,  0), epsilon), is(true));
		assertThat(new BigExp(6.66, -1).lessThan(new BigExp(6.66, +1), epsilon), is(true));
		assertThat(new BigExp(6.66,  0).lessThan(new BigExp(6.66, -1), epsilon), is(false));
		assertThat(new BigExp(6.66,  0).lessThan(new BigExp(6.66,  0), epsilon), is(false));
		assertThat(new BigExp(6.66,  0).lessThan(new BigExp(6.66, +1), epsilon), is(true));
		assertThat(new BigExp(6.66, +1).lessThan(new BigExp(6.66, -1), epsilon), is(false));
		assertThat(new BigExp(6.66, +1).lessThan(new BigExp(6.66,  0), epsilon), is(false));
		assertThat(new BigExp(6.66, +1).lessThan(new BigExp(6.66, +1), epsilon), is(false));

		// seen in the wild, caused by setting exponent of infinity to zero when normalizing
		assertThat(new BigExp(6.66, +1).lessThan(new BigExp(Double.POSITIVE_INFINITY)), is (true));
		assertThat(new BigExp(5.852204485365044,7).lessThan(new BigExp(Double.POSITIVE_INFINITY)), is(true));
	}

	@Test
	public void log() {

		// check accuracy to 16 significant digits (that's about all double can represent)
		var mc = new MathContext(16, RoundingMode.HALF_UP);
		BiFunction<Double,Integer,Double> be = (fp, exp) -> new BigExp(fp, exp).log();
		Function<String,Double> bd = val -> BigDecimalMath.log10(new BigDecimal(val), mc).doubleValue();

		assertThat(be.apply(1.2345, 10), isAbsolutely(bd.apply("1.2345e10"), 1e-12));
		assertThat(be.apply(1.2345, -10), isAbsolutely(bd.apply("1.2345e-10"), 1e-12));
	}

	@Test
	public void ln() {

		// check accuracy to 16 significant digits (that's about all double can represent)
		var mc = new MathContext(16, RoundingMode.HALF_UP);
		BiFunction<Double,Integer,Double> be = (fp, exp) -> new BigExp(fp, exp).ln();
		Function<String,Double> bd = val -> BigDecimalMath.log(new BigDecimal(val), mc).doubleValue();

		assertThat(be.apply(1.2345, 10), isAbsolutely(bd.apply("1.2345e10"), 1e-12));
		assertThat(be.apply(1.2345, -10), isAbsolutely(bd.apply("1.2345e-10"), 1e-12));
	}

	@Test
	public void pow10() {

		// check accuracy to 16 significant digits (that's about all double can represent)
		var mc = new MathContext(16, RoundingMode.HALF_UP);
		Function<Double,String> be = val -> String.format("%.12e", BigExp.pow10(val).toBigDecimal(mc));
		Function<String,String> bd = val -> String.format("%.12e", BigDecimalMath.pow(new BigDecimal("10"), new BigDecimal(val), mc));

		assertThat(be.apply(1.2345), is(bd.apply("1.2345")));
		assertThat(be.apply(12.345), is(bd.apply("12.345")));
		assertThat(be.apply(123.45), is(bd.apply("123.45")));
		assertThat(be.apply(308.0), is(bd.apply("308.0")));
		assertThat(be.apply(309.0), is(bd.apply("309.0")));
		assertThat(be.apply(1234.5), is(bd.apply("1234.5")));
		assertThat(be.apply(12345.0), is(bd.apply("12345")));

		assertThat(be.apply(-1.2345), is(bd.apply("-1.2345")));
		assertThat(be.apply(-12.345), is(bd.apply("-12.345")));
		assertThat(be.apply(-123.45), is(bd.apply("-123.45")));
		assertThat(be.apply(-300.0), is(bd.apply("-300.0")));
		assertThat(be.apply(-323.0), is(bd.apply("-323.0")));
		assertThat(be.apply(-324.0), is(bd.apply("-324.0")));
		assertThat(be.apply(-325.0), is(bd.apply("-325.0")));
		assertThat(be.apply(-1234.5), is(bd.apply("-1234.5")));
		assertThat(be.apply(-12345.0), is(bd.apply("-12345.0")));
	}

	@Test
	public void exp() {

		// check accuracy to 10 significant digits (the huger magnitude numbers are a little less precise)
		var mc = new MathContext(16, RoundingMode.HALF_UP);
		Function<Double,String> be = val -> String.format("%.10e", BigExp.exp(val).toBigDecimal(mc));
		Function<String,String> bd = val -> String.format("%.10e", BigDecimalMath.exp(new BigDecimal(val), mc));

		assertThat(be.apply(1.2345), is(bd.apply("1.2345")));
		assertThat(be.apply(12.345), is(bd.apply("12.345")));
		assertThat(be.apply(123.45), is(bd.apply("123.45")));
		assertThat(be.apply(308.0), is(bd.apply("308.0")));
		assertThat(be.apply(309.0), is(bd.apply("309.0")));
		assertThat(be.apply(1234.5), is(bd.apply("1234.5")));
		assertThat(be.apply(12345.0), is(bd.apply("12345.0")));

		assertThat(be.apply(-1.2345), is(bd.apply("-1.2345")));
		assertThat(be.apply(-12.345), is(bd.apply("-12.345")));
		assertThat(be.apply(-123.45), is(bd.apply("-123.45")));
		assertThat(be.apply(-300.0), is(bd.apply("-300.0")));
		assertThat(be.apply(-323.0), is(bd.apply("-323.0")));
		assertThat(be.apply(-324.0), is(bd.apply("-324.0")));
		assertThat(be.apply(-325.0), is(bd.apply("-325.0")));
		assertThat(be.apply(-1234.5), is(bd.apply("-1234.5")));
		assertThat(be.apply(-12345.0), is(bd.apply("-12345.0")));
	}

	private interface BigExpPow {
		String apply(double fp, int exp, double pow);
	}

	@Test
	public void pow() {

		// check accuracy to 16 significant digits (that's about all double can represent)
		var mc = new MathContext(16, RoundingMode.HALF_UP);
		BigExpPow be = (fp, exp, pow) -> {
			var val = new BigExp(fp, exp);
			val.pow(pow);
			return String.format("%.12e", val.toBigDecimal(mc));
		};
		BiFunction<String,Double,String> bd = (val, pow) -> String.format("%.12e", BigDecimalMath.pow(new BigDecimal(val), new BigDecimal(pow), mc));

		assertThat(be.apply(1.0, 0, 1.0), is(bd.apply("1.0", 1.0)));

		assertThat(be.apply(1.0, 0, 0.0), is(bd.apply("1.0", 0.0)));
		assertThat(be.apply(1.2345, 34, 0.0), is(bd.apply("1.2345e34", 0.0)));

		assertThat(be.apply(1.2345, 34, 3.0), is(bd.apply("1.2345e34", 3.0)));
		assertThat(be.apply(1.2345, 34, 3.7), is(bd.apply("1.2345e34", 3.7)));
		assertThat(be.apply(1.2345, 34, -3.7), is(bd.apply("1.2345e34", -3.7)));

		assertThat(be.apply(1.2345, 34, 0.456), is(bd.apply("1.2345e34", 0.456)));
		assertThat(be.apply(1.2345, 34, -0.456), is(bd.apply("1.2345e34", -0.456)));

		// huge numbers work too, but there's much less precision
		//assertThat(be.apply(1.2345, 6543, 204.6), is(bd.apply("1.2345e6543", 204.6)));
	}
}
