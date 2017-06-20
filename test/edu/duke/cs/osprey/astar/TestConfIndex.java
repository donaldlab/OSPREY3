package edu.duke.cs.osprey.astar;

import static edu.duke.cs.osprey.astar.Matchers.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

import edu.duke.cs.osprey.astar.conf.ConfIndex;

public class TestConfIndex {
	
	@Test
	public void isDefinedUndefinedRoot() {
		
		ConfIndex index = makeRoot5();
		index.numDefined = 0;
		index.numUndefined = 5;
		set(index.undefinedPos, 0, 1, 2, 3, 4);
		
		assertThat(index.isDefined(0), is(false));
		assertThat(index.isDefined(1), is(false));
		assertThat(index.isDefined(2), is(false));
		assertThat(index.isDefined(3), is(false));
		assertThat(index.isDefined(4), is(false));
		
		assertThat(index.isUndefined(0), is(true));
		assertThat(index.isUndefined(1), is(true));
		assertThat(index.isUndefined(2), is(true));
		assertThat(index.isUndefined(3), is(true));
		assertThat(index.isUndefined(4), is(true));
	}
	
	@Test
	public void isDefinedUndefinedChild() {
		
		ConfIndex index = new ConfIndex(5);
		index.numDefined = 3;
		set(index.definedPos, 0, 3, 4);
		set(index.definedRCs, 1, 2, 3);
		index.numUndefined = 2;
		set(index.undefinedPos, 1, 2);
		
		assertThat(index.isDefined(0), is(true));
		assertThat(index.isDefined(1), is(false));
		assertThat(index.isDefined(2), is(false));
		assertThat(index.isDefined(3), is(true));
		assertThat(index.isDefined(4), is(true));
		
		assertThat(index.isUndefined(0), is(false));
		assertThat(index.isUndefined(1), is(true));
		assertThat(index.isUndefined(2), is(true));
		assertThat(index.isUndefined(3), is(false));
		assertThat(index.isUndefined(4), is(false));
	}
	
	@Test
	public void isDefinedUndefinedWithGarbage() {
		
		ConfIndex index = new ConfIndex(5);
		index.numDefined = 3;
		set(index.definedPos, 0, 3, 4, 0, 1);
		set(index.definedRCs, 1, 2, 3, -1, -1);
		index.numUndefined = 2;
		set(index.undefinedPos, 1, 2, 4, 1, 3);
		
		assertThat(index.isDefined(0), is(true));
		assertThat(index.isDefined(1), is(false));
		assertThat(index.isDefined(2), is(false));
		assertThat(index.isDefined(3), is(true));
		assertThat(index.isDefined(4), is(true));
		
		assertThat(index.isUndefined(0), is(false));
		assertThat(index.isUndefined(1), is(true));
		assertThat(index.isUndefined(2), is(true));
		assertThat(index.isUndefined(3), is(false));
		assertThat(index.isUndefined(4), is(false));
	}
	
	@Test
	public void assignRoot() {
		
		// I AM ROOT
		ConfIndex index = makeRoot5()
			.assign(3, 5);
		
		assertThat(index.numPos, is(5));
		assertThat(index.numDefined, is(1));
		assertThat(index.definedPos, startsWith(3));
		assertThat(index.definedRCs, startsWith(5));
	}
	
	@Test
	public void assignChildBefore() {
		
		ConfIndex index = makeRoot5()
			.assign(3, 5)
			.assign(1, 6);
		
		assertThat(index.numPos, is(5));
		assertThat(index.numDefined, is(2));
		assertThat(index.definedPos, startsWith(1, 3));
		assertThat(index.definedRCs, startsWith(6, 5));
	}
	
	@Test
	public void assignChildAfter() {
		
		ConfIndex index = makeRoot5()
			.assign(3, 5)
			.assign(4, 6);
		
		assertThat(index.numPos, is(5));
		assertThat(index.numDefined, is(2));
		assertThat(index.definedPos, startsWith(3, 4));
		assertThat(index.definedRCs, startsWith(5, 6));
	}
	
	@Test
	public void assignAll01234() {
		
		ConfIndex index = makeRoot5()
			.assign(0, 4)
			.assign(1, 5)
			.assign(2, 6)
			.assign(3, 7)
			.assign(4, 8);
		
		assertThat(index.numPos, is(5));
		assertThat(index.numDefined, is(5));
		assertThat(index.definedPos, startsWith(0, 1, 2, 3, 4));
		assertThat(index.definedRCs, startsWith(4, 5, 6, 7, 8));
	}
	
	@Test
	public void assignAll43210() {
		
		ConfIndex index = makeRoot5()
			.assign(4, 8)
			.assign(3, 7)
			.assign(2, 6)
			.assign(1, 5)
			.assign(0, 4);
		
		assertThat(index.numPos, is(5));
		assertThat(index.numDefined, is(5));
		assertThat(index.definedPos, startsWith(0, 1, 2, 3, 4));
		assertThat(index.definedRCs, startsWith(4, 5, 6, 7, 8));
	}
	
	private ConfIndex makeRoot5() {
		ConfIndex confIndex = new ConfIndex(5);
		confIndex.numDefined = 0;
		confIndex.numUndefined = 5;
		set(confIndex.undefinedPos, 0, 1, 2, 3, 4);
		return confIndex;
	}
	
	private void set(int[] dest, int ... src) {
		for (int i=0; i<src.length; i++) {
			dest[i] = src[i];
		}
	}
}
