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

package edu.duke.cs.osprey.astar;

import java.io.Serializable;
import java.util.Arrays;

public class SlimAStarNode implements AStarNode, Serializable {
	
	private static final long serialVersionUID = -8612409900432350474L;

	public static class Factory implements AStarNode.Factory<SlimAStarNode>, Serializable {
		
		private static final long serialVersionUID = 4759571366344935185L;
		private int numPos;
		
		public Factory(int numPos) {
			this.numPos = numPos;
		}
		
		@Override
		public SlimAStarNode makeRoot() {
			return new SlimAStarNode(numPos);
		}

		@Override
		public SlimAStarNode make(SlimAStarNode parent, int pos, int rc) {
			return new SlimAStarNode(numPos, parent, pos, rc);
		}
	}
	
	public static class Link implements Comparable<Link>, Serializable {
		
		private static final long serialVersionUID = -6371596109692467264L;
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
	private short numPos; // TODO: any way to get rid of this?
	private short level; // NOTE: could get rid of this and have getLevel() compute it
	private Link link;
	private double gscore;
	private double hscore;
	
	private SlimAStarNode(int numPos, int level, Link link) {
		assert (numPos <= Short.MAX_VALUE);
		assert (level <= Short.MAX_VALUE);
		this.numPos = (short)numPos;
		this.level = (short)level;
		this.link = link;
		gscore = Double.NaN;
		hscore = Double.NaN;
	}
	
	private SlimAStarNode(int numPos) {
		this(numPos, 0, new Link());
	}

	private SlimAStarNode(int numPos, SlimAStarNode parent, int assignedPos, int assignedRc) {
		this(numPos, parent.level + 1, new Link(parent.link, assignedPos, assignedRc));
	}
	
	public Link getLink() {
		return link;
	}
	
	@Override
	public int compareTo(AStarNode other) {
		return Double.compare(getScore(), other.getScore());
	}

	@Override
	public int[] getNodeAssignments() {
		
		// NOTE: don't call this too much, it's very slow
		// try to use PairwiseConfTree.splitPositions() where you can
		
		int[] conf = new int[numPos];
		Arrays.fill(conf, -1);
		Link link = this.link;
		while (!link.isRoot()) {
			conf[link.getPos()] = link.getRC();
			link = link.getParent();
		}
		return conf;
	}
	
	@Override
	public double getScore() {
		return gscore + hscore;
	}
	
	public void setScore(double val) {
		throw new UnsupportedOperationException("use setGScore() and setHScore() instead");
	}

	public double getGScore() {
		return gscore;
	}
	public void setGScore(double val) {
		gscore = val;
	}
	
	public double getHScore() {
		return hscore;
	}
	public void setHScore(double val) {
		hscore = val;
	}

	@Override
	public boolean scoreNeedsRefinement() {
		return false;
	}

	@Override
	public void setScoreNeedsRefinement(boolean val) {
		// do nothing
	}
	
	@Override
	public int getLevel() {
		return level;
	}

	@Override
	public boolean isFullyDefined() {
		return getLevel() == numPos;
	}
}
