package edu.duke.cs.osprey.tools;

import static org.junit.Assert.*;

import org.junit.Test;

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
}
