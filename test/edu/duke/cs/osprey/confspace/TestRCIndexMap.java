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
