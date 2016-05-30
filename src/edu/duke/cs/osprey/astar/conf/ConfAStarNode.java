package edu.duke.cs.osprey.astar.conf;

import java.util.Arrays;

public class ConfAStarNode implements Comparable<ConfAStarNode> {

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
	
	private ConfAStarNode(int level, Link link) {
		assert (level <= Short.MAX_VALUE);
		this.level = (short)level;
		this.link = link;
		gscore = Double.NaN;
		hscore = Double.NaN;
	}
	
	public ConfAStarNode() {
		this(0, new Link());
	}

	public ConfAStarNode(ConfAStarNode parent, int assignedPos, int assignedRc) {
		this(parent.level + 1, new Link(parent.link, assignedPos, assignedRc));
	}
	
	public Link getLink() {
		return link;
	}
	
	@Override
	public int compareTo(ConfAStarNode other) {
		return Double.compare(getScore(), other.getScore());
	}

	public void getConf(int[] conf) {
		Arrays.fill(conf, -1);
		Link link = this.link;
		while (!link.isRoot()) {
			conf[link.getPos()] = link.getRC();
			link = link.getParent();
		}
	}
	
	public double getScore() {
		return gscore + hscore;
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

	public int getLevel() {
		return level;
	}
}
