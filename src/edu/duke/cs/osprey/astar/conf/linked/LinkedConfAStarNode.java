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

package edu.duke.cs.osprey.astar.conf.linked;

import java.util.Arrays;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;

public class LinkedConfAStarNode implements ConfAStarNode {
	
	public static class Link implements Comparable<Link> {
		
		// NOTE: try to keep storage here as small as possible
		// we expect to have millions of nodes in memory
		private Link parent;
		private short pos;
		private short rc;
		
		public Link() {
			this(null, -1, -1);
		}
		
		public Link(Link parent, int pos, int rc) {
			assert (pos <= Short.MAX_VALUE);
			assert (rc <= Short.MAX_VALUE);
			this.parent = parent;
			this.pos = (short)pos;
			this.rc = (short)rc;
		}
		
		public Link getParent() {
			return parent;
		}
		
		public int getPos() {
			return pos;
		}
		
		public int getRC() {
			return rc;
		}
		
		public boolean isRoot() {
			return parent == null;
		}

		@Override
		public int compareTo(Link other) {
			return pos - other.pos;
		}
	}
	
	// NOTE: try to keep storage here as small as possible
	// we expect to have millions of nodes in memory
	private short level; // NOTE: could get rid of this and have getLevel() compute it
	private Link link;
	private double gscore;
	private double hscore;
	
	private LinkedConfAStarNode(int level, Link link) {
		assert (level <= Short.MAX_VALUE);
		this.level = (short)level;
		this.link = link;
		gscore = Double.NaN;
		hscore = Double.NaN;
	}
	
	public LinkedConfAStarNode() {
		this(0, new Link());
	}
	
	@Override
	public LinkedConfAStarNode assign(int pos, int rc) {
		return new LinkedConfAStarNode(level + 1, new Link(link, pos, rc));
	}
	
	public Link getLink() {
		return link;
	}
	
	@Override
	public void getConf(int[] conf) {
		Arrays.fill(conf, -1);
		Link link = this.link;
		while (!link.isRoot()) {
			conf[link.getPos()] = link.getRC();
			link = link.getParent();
		}
	}
	
	@Override
	public double getGScore() {
		return gscore;
	}
	
	@Override
	public void setGScore(double val) {
		gscore = val;
	}
	
	@Override
	public double getHScore() {
		return hscore;
	}
	
	@Override
	public void setHScore(double val) {
		hscore = val;
	}

	@Override
	public int getLevel() {
		return level;
	}
	
	@Override
	public void index(ConfIndex index) {
		
		// is this node already indexed?
		if (index.node == this) {
			return;
		}
		index.node = this;
		
		// use local vars so the (JIT)compiler can use stack/registers instead of field accesses
		int numPos = index.numPos;
		int numDefined = 0;
		int[] dpos = index.definedPos;
		int[] rcs = index.definedRCs;
		int numUndefined = 0;
		int[] upos = index.undefinedPos;
		
		// split conformation into defined and undefined positions
		
		// do one pass through the link chain to get the defined positions
		LinkedConfAStarNode.Link link = this.link;
		while (!(link.isRoot())) {
			dpos[numDefined] = link.getPos();
			rcs[numDefined] = link.getRC();
			numDefined++;
			link = link.getParent();
		}
		
		// sort the defined positions using a simple insertion sort
		// assignments arrays are always small (n << 100), so insertion sort should be fast enough
		// NOTE: we need to sort two arrays simultaneously, so we can't use any library sorts
		{
			for (int i=1; i<numDefined; i++) {
				
				int tempPos = dpos[i];
				int tempRC = rcs[i];
				
				int j;
				for (j=i; j>=1 && tempPos < dpos[j-1]; j--) {
					dpos[j] = dpos[j-1];
					rcs[j] = rcs[j-1];
				}
				dpos[j] = tempPos;
				rcs[j] = tempRC;
			}
		}
		
		// now figure out the undefined positions
		if (numDefined == 0) {
			
			// all undefined
			numUndefined = numPos;
			for (int pos=0; pos<numPos; pos++) {
				upos[pos] = pos;
			}
			
		} else {
			
			int i = 0;
			for (int pos=0; pos<numPos; pos++) {
				
				// does this pos match the next defined pos?
				if (i < numDefined && pos == dpos[i]) {
					i++;
				} else {
					upos[numUndefined] = pos;
					numUndefined++;
				}
			}
		}
		
		assert (numDefined + numUndefined == numPos);
		
		// copy vars back to the index
		index.numDefined = numDefined;
		index.numUndefined = numUndefined;
	}
}
