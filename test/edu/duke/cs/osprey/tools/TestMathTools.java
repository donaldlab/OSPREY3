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

import static org.junit.Assert.*;

import org.junit.Test;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.Iterator;

import static org.hamcrest.Matchers.*;

public class TestMathTools {
	
	@Test
	public void divUp() {
		
		assertThat(MathTools.divUp(0, 1), is(0));
		assertThat(MathTools.divUp(0, 2), is(0));
		assertThat(MathTools.divUp(0, 3), is(0));
		assertThat(MathTools.divUp(0, 4), is(0));
		
		assertThat(MathTools.divUp(1, 1), is(1));
		assertThat(MathTools.divUp(1, 2), is(1));
		assertThat(MathTools.divUp(1, 3), is(1));
		assertThat(MathTools.divUp(1, 4), is(1));
		
		assertThat(MathTools.divUp(2, 1), is(2));
		assertThat(MathTools.divUp(2, 2), is(1));
		assertThat(MathTools.divUp(2, 3), is(1));
		assertThat(MathTools.divUp(2, 4), is(1));
		
		assertThat(MathTools.divUp(3, 1), is(3));
		assertThat(MathTools.divUp(3, 2), is(2));
		assertThat(MathTools.divUp(3, 3), is(1));
		assertThat(MathTools.divUp(3, 4), is(1));
		
		assertThat(MathTools.divUp(4, 1), is(4));
		assertThat(MathTools.divUp(4, 2), is(2));
		assertThat(MathTools.divUp(4, 3), is(2));
		assertThat(MathTools.divUp(4, 4), is(1));
		
		assertThat(MathTools.divUp(5, 1), is(5));
		assertThat(MathTools.divUp(5, 2), is(3));
		assertThat(MathTools.divUp(5, 3), is(2));
		assertThat(MathTools.divUp(5, 4), is(2));
	}
	
	@Test
	public void roundUpMultiple() {
		
		assertThat(MathTools.roundUpToMultiple(0, 1), is(0));
		assertThat(MathTools.roundUpToMultiple(1, 1), is(1));
		assertThat(MathTools.roundUpToMultiple(2, 1), is(2));
		assertThat(MathTools.roundUpToMultiple(3, 1), is(3));
		
		assertThat(MathTools.roundUpToMultiple(0, 2), is(0));
		assertThat(MathTools.roundUpToMultiple(1, 2), is(2));
		assertThat(MathTools.roundUpToMultiple(2, 2), is(2));
		assertThat(MathTools.roundUpToMultiple(3, 2), is(4));
		assertThat(MathTools.roundUpToMultiple(4, 2), is(4));
		
		assertThat(MathTools.roundUpToMultiple(0, 4), is(0));
		assertThat(MathTools.roundUpToMultiple(1, 4), is(4));
		assertThat(MathTools.roundUpToMultiple(2, 4), is(4));
		assertThat(MathTools.roundUpToMultiple(3, 4), is(4));
		assertThat(MathTools.roundUpToMultiple(4, 4), is(4));
		assertThat(MathTools.roundUpToMultiple(5, 4), is(8));
		assertThat(MathTools.roundUpToMultiple(6, 4), is(8));
		assertThat(MathTools.roundUpToMultiple(7, 4), is(8));
		assertThat(MathTools.roundUpToMultiple(8, 4), is(8));
		assertThat(MathTools.roundUpToMultiple(9, 4), is(12));
	}

	@Test
	public void powerset() {

		// n = 1
		assertThat(MathTools.powerset(Arrays.asList(1)), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1)
		)));

		// n = 2
		assertThat(MathTools.powerset(Arrays.asList(1, 2)), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1),
			Arrays.asList(2),
			Arrays.asList(1, 2)
		)));

		// n = 3
		assertThat(MathTools.powerset(Arrays.asList(1, 2, 3)), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1),
			Arrays.asList(2),
			Arrays.asList(1, 2),
			Arrays.asList(3),
			Arrays.asList(1, 3),
			Arrays.asList(2, 3),
			Arrays.asList(1, 2, 3)
		)));

		// n = 4
		assertThat(MathTools.powerset(Arrays.asList(1, 2, 3, 4)), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1),
			Arrays.asList(2),
			Arrays.asList(1, 2),
			Arrays.asList(3),
			Arrays.asList(1, 3),
			Arrays.asList(2, 3),
			Arrays.asList(1, 2, 3),
			Arrays.asList(4),
			Arrays.asList(1, 4),
			Arrays.asList(2, 4),
			Arrays.asList(1, 2, 4),
			Arrays.asList(3, 4),
			Arrays.asList(1, 3, 4),
			Arrays.asList(2, 3, 4),
			Arrays.asList(1, 2, 3, 4)
		)));
	}

	@Test
	public void powersetUpTo() {

		// n = 1, size = 0
		assertThat(MathTools.powersetUpTo(Arrays.asList(1), 0), is(Arrays.asList(
			Arrays.asList()
		)));

		// n = 1, size = 1
		assertThat(MathTools.powersetUpTo(Arrays.asList(1), 1), is(MathTools.powerset(Arrays.asList(1))));

		// n = 2, size = 0
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2), 0), is(Arrays.asList(
			Arrays.asList()
		)));

		// n = 2, size = 1
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2), 1), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1),
			Arrays.asList(2)
		)));

		// n = 2, size = 2
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2), 2), is(MathTools.powerset(Arrays.asList(1, 2))));

		// n = 3, size = 0
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2, 3), 0), is(Arrays.asList(
			Arrays.asList()
		)));

		// n = 3, size = 1
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2, 3), 1), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1),
			Arrays.asList(2),
			Arrays.asList(3)
		)));

		// n = 3, size = 2
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2, 3), 2), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1),
			Arrays.asList(2),
			Arrays.asList(1, 2),
			Arrays.asList(3),
			Arrays.asList(1, 3),
			Arrays.asList(2, 3)
		)));

		// n = 3, size = 3
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2, 3), 3), is(MathTools.powerset(Arrays.asList(1, 2, 3))));

		// n = 4, size = 0
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2, 3, 4), 0), is(Arrays.asList(
			Arrays.asList()
		)));

		// n = 4, size = 1
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2, 3, 4), 1), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1),
			Arrays.asList(2),
			Arrays.asList(3),
			Arrays.asList(4)
		)));

		// n = 4, size = 2
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2, 3, 4), 2), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1),
			Arrays.asList(2),
			Arrays.asList(1, 2),
			Arrays.asList(3),
			Arrays.asList(1, 3),
			Arrays.asList(2, 3),
			Arrays.asList(4),
			Arrays.asList(1, 4),
			Arrays.asList(2, 4),
			Arrays.asList(3, 4)
		)));

		// n = 4, size = 3
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2, 3, 4), 3), is(Arrays.asList(
			Arrays.asList(),
			Arrays.asList(1),
			Arrays.asList(2),
			Arrays.asList(1, 2),
			Arrays.asList(3),
			Arrays.asList(1, 3),
			Arrays.asList(2, 3),
			Arrays.asList(1, 2, 3),
			Arrays.asList(4),
			Arrays.asList(1, 4),
			Arrays.asList(2, 4),
			Arrays.asList(1, 2, 4),
			Arrays.asList(3, 4),
			Arrays.asList(1, 3, 4),
			Arrays.asList(2, 3, 4)
		)));

		// n = 4, size = 4
		assertThat(MathTools.powersetUpTo(Arrays.asList(1, 2, 3, 4), 4), is(MathTools.powerset(Arrays.asList(1, 2, 3, 4))));
	}

	@Test
	@SuppressWarnings("unchecked")
	public void cartesianProduct() {

		// empty
		assertThat(MathTools.cartesianProduct(Arrays.asList()).iterator().hasNext(), is(false));

		// 0
		assertThat(MathTools.cartesianProduct(Arrays.asList(
			Arrays.asList()
		)).iterator().hasNext(), is(false));

		// 2 x 0
		assertThat(MathTools.cartesianProduct(Arrays.asList(
			Arrays.asList("a", "b"),
			Arrays.asList()
		)).iterator().hasNext(), is(false));

		// 2 x 2
		assertThat(MathTools.cartesianProduct(Arrays.asList(
			Arrays.asList("a", "b"),
			Arrays.asList("c", "d")
		)), contains(
			Arrays.asList("a", "c"),
			Arrays.asList("a", "d"),
			Arrays.asList("b", "c"),
			Arrays.asList("b", "d")
		));

		// 2 x 3
		assertThat(MathTools.cartesianProduct(Arrays.asList(
			Arrays.asList("a", "b", "c"),
			Arrays.asList("d", "e")
		)), contains(
			Arrays.asList("a", "d"),
			Arrays.asList("a", "e"),
			Arrays.asList("b", "d"),
			Arrays.asList("b", "e"),
			Arrays.asList("c", "d"),
			Arrays.asList("c", "e")
		));
	}

	@Test
	public void isLessThan() {

		assertThat(MathTools.isLessThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(Double.POSITIVE_INFINITY)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(5.0)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(0.0)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(-5.0)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(false));

		assertThat(MathTools.isLessThan(MathTools.biggen(5.0), MathTools.biggen(Double.POSITIVE_INFINITY)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(5.0), MathTools.biggen(5.0)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(5.0), MathTools.biggen(0.0)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(5.0), MathTools.biggen(-5.0)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(5.0), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(false));

		assertThat(MathTools.isLessThan(MathTools.biggen(0.0), MathTools.biggen(Double.POSITIVE_INFINITY)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(0.0), MathTools.biggen(5.0)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(0.0), MathTools.biggen(0.0)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(0.0), MathTools.biggen(-5.0)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(0.0), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(false));

		assertThat(MathTools.isLessThan(MathTools.biggen(-5.0), MathTools.biggen(Double.POSITIVE_INFINITY)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(-5.0), MathTools.biggen(5.0)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(-5.0), MathTools.biggen(0.0)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(-5.0), MathTools.biggen(-5.0)), is(false));
		assertThat(MathTools.isLessThan(MathTools.biggen(-5.0), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(false));

		assertThat(MathTools.isLessThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(Double.POSITIVE_INFINITY)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(5.0)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(0.0)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(-5.0)), is(true));
		assertThat(MathTools.isLessThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(false));
	}

	@Test
	public void isGreaterThan() {

		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(Double.POSITIVE_INFINITY)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(5.0)), is(true));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(0.0)), is(true));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(-5.0)), is(true));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.POSITIVE_INFINITY), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(true));

		assertThat(MathTools.isGreaterThan(MathTools.biggen(5.0), MathTools.biggen(Double.POSITIVE_INFINITY)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(5.0), MathTools.biggen(5.0)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(5.0), MathTools.biggen(0.0)), is(true));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(5.0), MathTools.biggen(-5.0)), is(true));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(5.0), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(true));

		assertThat(MathTools.isGreaterThan(MathTools.biggen(0.0), MathTools.biggen(Double.POSITIVE_INFINITY)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(0.0), MathTools.biggen(5.0)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(0.0), MathTools.biggen(0.0)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(0.0), MathTools.biggen(-5.0)), is(true));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(0.0), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(true));

		assertThat(MathTools.isGreaterThan(MathTools.biggen(-5.0), MathTools.biggen(Double.POSITIVE_INFINITY)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(-5.0), MathTools.biggen(5.0)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(-5.0), MathTools.biggen(0.0)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(-5.0), MathTools.biggen(-5.0)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(-5.0), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(true));

		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(Double.POSITIVE_INFINITY)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(5.0)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(0.0)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(-5.0)), is(false));
		assertThat(MathTools.isGreaterThan(MathTools.biggen(Double.NEGATIVE_INFINITY), MathTools.biggen(Double.NEGATIVE_INFINITY)), is(false));
	}

	public static BigDecimal bigAdd(double a, double b) {
		return MathTools.bigAdd(
			MathTools.biggen(a),
			MathTools.biggen(b),
			new MathContext(64, RoundingMode.HALF_UP)
		);
	}

	@Test
	public void bigAdd() {

		assertThat(bigAdd(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY), is(MathTools.BigPositiveInfinity));
		assertThat(bigAdd(Double.POSITIVE_INFINITY, 5.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigAdd(Double.POSITIVE_INFINITY, 0.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigAdd(Double.POSITIVE_INFINITY, -5.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigAdd(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY), is(MathTools.BigNaN));

		assertThat(bigAdd(5.0, Double.POSITIVE_INFINITY), is(MathTools.BigPositiveInfinity));
		assertThat(bigAdd(5.0, 5.0).doubleValue(), is(10.0));
		assertThat(bigAdd(5.0, 0.0).doubleValue(), is(5.0));
		assertThat(bigAdd(5.0, -5.0).doubleValue(), is(0.0));
		assertThat(bigAdd(5.0, Double.NEGATIVE_INFINITY), is(MathTools.BigNegativeInfinity));

		assertThat(bigAdd(0.0, Double.POSITIVE_INFINITY), is(MathTools.BigPositiveInfinity));
		assertThat(bigAdd(0.0, 5.0).doubleValue(), is(5.0));
		assertThat(bigAdd(0.0, 0.0).doubleValue(), is(0.0));
		assertThat(bigAdd(0.0, -5.0).doubleValue(), is(-5.0));
		assertThat(bigAdd(0.0, Double.NEGATIVE_INFINITY), is(MathTools.BigNegativeInfinity));

		assertThat(bigAdd(-5.0, Double.POSITIVE_INFINITY), is(MathTools.BigPositiveInfinity));
		assertThat(bigAdd(-5.0, 5.0).doubleValue(), is(0.0));
		assertThat(bigAdd(-5.0, 0.0).doubleValue(), is(-5.0));
		assertThat(bigAdd(-5.0, -5.0).doubleValue(), is(-10.0));
		assertThat(bigAdd(-5.0, Double.NEGATIVE_INFINITY), is(MathTools.BigNegativeInfinity));

		assertThat(bigAdd(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY), is(MathTools.BigNaN));
		assertThat(bigAdd(Double.NEGATIVE_INFINITY, 5.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigAdd(Double.NEGATIVE_INFINITY, 0.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigAdd(Double.NEGATIVE_INFINITY, -5.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigAdd(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY), is(MathTools.BigNegativeInfinity));
	}

	public static BigDecimal bigSubtract(double a, double b) {
		return MathTools.bigSubtract(
			MathTools.biggen(a),
			MathTools.biggen(b),
			new MathContext(64, RoundingMode.HALF_UP)
		);
	}

	@Test
	public void bigSubtract() {

		assertThat(bigSubtract(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY), is(MathTools.BigNaN));
		assertThat(bigSubtract(Double.POSITIVE_INFINITY, 5.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigSubtract(Double.POSITIVE_INFINITY, 0.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigSubtract(Double.POSITIVE_INFINITY, -5.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigSubtract(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY), is(MathTools.BigPositiveInfinity));

		assertThat(bigSubtract(5.0, Double.POSITIVE_INFINITY), is(MathTools.BigNegativeInfinity));
		assertThat(bigSubtract(5.0, 5.0).doubleValue(), is(0.0));
		assertThat(bigSubtract(5.0, 0.0).doubleValue(), is(5.0));
		assertThat(bigSubtract(5.0, -5.0).doubleValue(), is(10.0));
		assertThat(bigSubtract(5.0, Double.NEGATIVE_INFINITY), is(MathTools.BigPositiveInfinity));

		assertThat(bigSubtract(0.0, Double.POSITIVE_INFINITY), is(MathTools.BigNegativeInfinity));
		assertThat(bigSubtract(0.0, 5.0).doubleValue(), is(-5.0));
		assertThat(bigSubtract(0.0, 0.0).doubleValue(), is(0.0));
		assertThat(bigSubtract(0.0, -5.0).doubleValue(), is(5.0));
		assertThat(bigSubtract(0.0, Double.NEGATIVE_INFINITY), is(MathTools.BigPositiveInfinity));

		assertThat(bigSubtract(-5.0, Double.POSITIVE_INFINITY), is(MathTools.BigNegativeInfinity));
		assertThat(bigSubtract(-5.0, 5.0).doubleValue(), is(-10.0));
		assertThat(bigSubtract(-5.0, 0.0).doubleValue(), is(-5.0));
		assertThat(bigSubtract(-5.0, -5.0).doubleValue(), is(0.0));
		assertThat(bigSubtract(-5.0, Double.NEGATIVE_INFINITY), is(MathTools.BigPositiveInfinity));

		assertThat(bigSubtract(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY), is(MathTools.BigNegativeInfinity));
		assertThat(bigSubtract(Double.NEGATIVE_INFINITY, 5.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigSubtract(Double.NEGATIVE_INFINITY, 0.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigSubtract(Double.NEGATIVE_INFINITY, -5.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigSubtract(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY), is(MathTools.BigNaN));
	}

	public static BigDecimal bigMultiply(double a, double b) {
		return MathTools.bigMultiply(
			MathTools.biggen(a),
			MathTools.biggen(b),
			new MathContext(64, RoundingMode.HALF_UP)
		);
	}

	@Test
	public void bigMultiply() {

		assertThat(bigMultiply(Double.NaN, Double.NaN), is(MathTools.BigNaN));

		assertThat(bigMultiply(Double.POSITIVE_INFINITY, Double.NaN), is(MathTools.BigNaN));
		assertThat(bigMultiply(Double.NEGATIVE_INFINITY, Double.NaN), is(MathTools.BigNaN));
		assertThat(bigMultiply(-1.0, Double.NaN), is(MathTools.BigNaN));
		assertThat(bigMultiply(0.0, Double.NaN), is(MathTools.BigNaN));
		assertThat(bigMultiply(1.0, Double.NaN), is(MathTools.BigNaN));

		assertThat(bigMultiply(Double.NaN, Double.POSITIVE_INFINITY), is(MathTools.BigNaN));
		assertThat(bigMultiply(Double.NaN, Double.NEGATIVE_INFINITY), is(MathTools.BigNaN));
		assertThat(bigMultiply(Double.NaN, -1.0), is(MathTools.BigNaN));
		assertThat(bigMultiply(Double.NaN, 0.0), is(MathTools.BigNaN));
		assertThat(bigMultiply(Double.NaN, 1.0), is(MathTools.BigNaN));

		assertThat(bigMultiply(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY), is(MathTools.BigPositiveInfinity));
		assertThat(bigMultiply(Double.POSITIVE_INFINITY, 5.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigMultiply(Double.POSITIVE_INFINITY, 0.0).doubleValue(), is(0.0));
		assertThat(bigMultiply(Double.POSITIVE_INFINITY, -5.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigMultiply(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY), is(MathTools.BigNaN));

		assertThat(bigMultiply(5.0, Double.POSITIVE_INFINITY), is(MathTools.BigPositiveInfinity));
		assertThat(bigMultiply(5.0, 5.0).doubleValue(), is(25.0));
		assertThat(bigMultiply(5.0, 0.0).doubleValue(), is(0.0));
		assertThat(bigMultiply(5.0, -5.0).doubleValue(), is(-25.0));
		assertThat(bigMultiply(5.0, Double.NEGATIVE_INFINITY), is(MathTools.BigNegativeInfinity));

		assertThat(bigMultiply(0.0, Double.POSITIVE_INFINITY).doubleValue(), is(0.0));
		assertThat(bigMultiply(0.0, 5.0).doubleValue(), is(0.0));
		assertThat(bigMultiply(0.0, 0.0).doubleValue(), is(0.0));
		assertThat(bigMultiply(0.0, -5.0).doubleValue(), is(0.0));
		assertThat(bigMultiply(0.0, Double.NEGATIVE_INFINITY).doubleValue(), is(0.0));

		assertThat(bigMultiply(-5.0, Double.POSITIVE_INFINITY), is(MathTools.BigNegativeInfinity));
		assertThat(bigMultiply(-5.0, 5.0).doubleValue(), is(-25.0));
		assertThat(bigMultiply(-5.0, 0.0).doubleValue(), is(0.0));
		assertThat(bigMultiply(-5.0, -5.0).doubleValue(), is(25.0));
		assertThat(bigMultiply(-5.0, Double.NEGATIVE_INFINITY), is(MathTools.BigPositiveInfinity));

		assertThat(bigMultiply(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY), is(MathTools.BigNaN));
		assertThat(bigMultiply(Double.NEGATIVE_INFINITY, 5.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigMultiply(Double.NEGATIVE_INFINITY, 0.0).doubleValue(), is(0.0));
		assertThat(bigMultiply(Double.NEGATIVE_INFINITY, -5.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigMultiply(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY), is(MathTools.BigPositiveInfinity));
	}

	public static BigDecimal bigDivide(double a, double b) {
		return MathTools.bigDivide(
			MathTools.biggen(a),
			MathTools.biggen(b),
			new MathContext(64, RoundingMode.HALF_UP)
		);
	}

	@Test
	public void bigDivide() {

		assertThat(bigDivide(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY), is(MathTools.BigNaN));
		assertThat(bigDivide(Double.POSITIVE_INFINITY, 5.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigDivide(Double.POSITIVE_INFINITY, 0.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigDivide(Double.POSITIVE_INFINITY, -5.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigDivide(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY), is(MathTools.BigNaN));

		assertThat(bigDivide(5.0, Double.POSITIVE_INFINITY).doubleValue(), is(0.0));
		assertThat(bigDivide(5.0, 5.0).doubleValue(), is(1.0));
		assertThat(bigDivide(5.0, 0.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigDivide(5.0, -5.0).doubleValue(), is(-1.0));
		assertThat(bigDivide(5.0, Double.NEGATIVE_INFINITY).doubleValue(), is(0.0));

		assertThat(bigDivide(0.0, Double.POSITIVE_INFINITY).doubleValue(), is(0.0));
		assertThat(bigDivide(0.0, 5.0).doubleValue(), is(0.0));
		assertThat(bigDivide(0.0, 0.0), is(MathTools.BigNaN));
		assertThat(bigDivide(0.0, -5.0).doubleValue(), is(0.0));
		assertThat(bigDivide(0.0, Double.NEGATIVE_INFINITY).doubleValue(), is(0.0));

		assertThat(bigDivide(-5.0, Double.POSITIVE_INFINITY).doubleValue(), is(0.0));
		assertThat(bigDivide(-5.0, 5.0).doubleValue(), is(-1.0));
		assertThat(bigDivide(-5.0, 0.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigDivide(-5.0, -5.0).doubleValue(), is(1.0));
		assertThat(bigDivide(-5.0, Double.NEGATIVE_INFINITY).doubleValue(), is(0.0));

		assertThat(bigDivide(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY), is(MathTools.BigNaN));
		assertThat(bigDivide(Double.NEGATIVE_INFINITY, 5.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigDivide(Double.NEGATIVE_INFINITY, 0.0), is(MathTools.BigNegativeInfinity));
		assertThat(bigDivide(Double.NEGATIVE_INFINITY, -5.0), is(MathTools.BigPositiveInfinity));
		assertThat(bigDivide(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY), is(MathTools.BigNaN));
	}

	public static double bigNegate(double a) {
		return MathTools.bigNegate(MathTools.biggen(a)).doubleValue();
	}

	@Test
	public void bigNegate() {
		assertThat(bigNegate(-1.0), is(1.0));
		assertThat(bigNegate(0.0), is(0.0));
		assertThat(bigNegate(1.0), is(-1.0));
		assertThat(bigNegate(Double.NaN), is(Double.NaN));
		assertThat(bigNegate(Double.NEGATIVE_INFINITY), is(Double.POSITIVE_INFINITY));
		assertThat(bigNegate(Double.POSITIVE_INFINITY), is(Double.NEGATIVE_INFINITY));
	}

	@Test
	public void gridIterable() {

		Iterator<int[]> iter = new MathTools.GridIterable(new int[] { 1, 2, 3 }).iterator();

		assertThat(iter.hasNext(), is(true));
		assertThat(iter.next(), is(new int[] { 0, 0, 0 }));

		assertThat(iter.hasNext(), is(true));
		assertThat(iter.next(), is(new int[] { 0, 1, 0 }));

		assertThat(iter.hasNext(), is(true));
		assertThat(iter.next(), is(new int[] { 0, 0, 1 }));

		assertThat(iter.hasNext(), is(true));
		assertThat(iter.next(), is(new int[] { 0, 1, 1 }));

		assertThat(iter.hasNext(), is(true));
		assertThat(iter.next(), is(new int[] { 0, 0, 2 }));

		assertThat(iter.hasNext(), is(true));
		assertThat(iter.next(), is(new int[] { 0, 1, 2 }));

		assertThat(iter.hasNext(), is(false));
	}
}
