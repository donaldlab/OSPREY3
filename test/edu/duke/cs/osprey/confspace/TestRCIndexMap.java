package edu.duke.cs.osprey.confspace;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

public class TestRCIndexMap {
	
	@Test
	public void testBase() {
		RCIndexMap map = new RCIndexMap(5);
		
		// old->new
		assertThat(map.oldToNew(0), contains(0));
		assertThat(map.oldToNew(1), contains(1));
		assertThat(map.oldToNew(2), contains(2));
		assertThat(map.oldToNew(3), contains(3));
		assertThat(map.oldToNew(4), contains(4));
		
		// new->old
		assertThat(map.newToOld(0), is(0));
		assertThat(map.newToOld(1), is(1));
		assertThat(map.newToOld(2), is(2));
		assertThat(map.newToOld(3), is(3));
		assertThat(map.newToOld(4), is(4));
	}
	
	@Test
	public void testRemove() {
		RCIndexMap map = new RCIndexMap(5);
		
		map.remove(2);
		
		// old->new
		assertThat(map.oldToNew(0), contains(0));
		assertThat(map.oldToNew(1), contains(1));
		assertThat(map.oldToNew(2).isEmpty(), is(true));
		assertThat(map.oldToNew(3), contains(3));
		assertThat(map.oldToNew(4), contains(4));
		
		// new->old
		assertThat(map.newToOld(0), is(0));
		assertThat(map.newToOld(1), is(1));
		assertThat(map.newToOld(2), is(nullValue()));
		assertThat(map.newToOld(3), is(3));
		assertThat(map.newToOld(4), is(4));
	}
	
	@Test
	public void testAdd() {
		RCIndexMap map = new RCIndexMap(5);
		
		map.add(2, 5);
		
		// old->new
		assertThat(map.oldToNew(0), contains(0));
		assertThat(map.oldToNew(1), contains(1));
		assertThat(map.oldToNew(2), contains(2, 5));
		assertThat(map.oldToNew(3), contains(3));
		assertThat(map.oldToNew(4), contains(4));
		
		// new->old
		assertThat(map.newToOld(0), is(0));
		assertThat(map.newToOld(1), is(1));
		assertThat(map.newToOld(2), is(2));
		assertThat(map.newToOld(3), is(3));
		assertThat(map.newToOld(4), is(4));
		assertThat(map.newToOld(5), is(nullValue()));
	}
	
	@Test
	public void testRemoveAdd() {
		RCIndexMap map = new RCIndexMap(5);
		
		map.remove(2);
		map.add(2, 2);
		
		// old->new
		assertThat(map.oldToNew(0), contains(0));
		assertThat(map.oldToNew(1), contains(1));
		assertThat(map.oldToNew(2), contains(2));
		assertThat(map.oldToNew(3), contains(3));
		assertThat(map.oldToNew(4), contains(4));
		
		// new->old
		assertThat(map.newToOld(0), is(0));
		assertThat(map.newToOld(1), is(1));
		assertThat(map.newToOld(2), is(nullValue()));
		assertThat(map.newToOld(3), is(3));
		assertThat(map.newToOld(4), is(4));
	}
}
