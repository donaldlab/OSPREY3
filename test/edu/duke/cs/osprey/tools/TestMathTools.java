package edu.duke.cs.osprey.tools;

import static org.junit.Assert.*;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

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
}
