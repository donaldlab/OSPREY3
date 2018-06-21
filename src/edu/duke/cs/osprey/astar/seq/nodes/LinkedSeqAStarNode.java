package edu.duke.cs.osprey.astar.seq.nodes;

import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.Sequence;


public class LinkedSeqAStarNode implements SeqAStarNode {

	public static class Link implements Comparable<Link> {

		// NOTE: try to keep storage here as small as possible
		// we expect to have millions of nodes in memory
		public final Link parent;
		public final short pos;
		public final short rt;

		public Link() {
			this(null, -1, -1);
		}

		public Link(Link parent, int pos, int rt) {
			assert (pos <= Short.MAX_VALUE);
			assert (rt <= Short.MAX_VALUE);
			this.parent = parent;
			this.pos = (short) pos;
			this.rt = (short)rt;
		}

		public boolean isRoot() {
			return parent == null;
		}

		@Override
		public int compareTo(Link other) {
			return this.pos - other.pos;
		}
	}

	// NOTE: try to keep storage here as small as possible
	// we expect to have millions of nodes in memory
	private short level; // NOTE: could get rid of this and have getLevel() compute it
	private Link link;
	private double gscore;
	private double hscore;
	private Object data;

	private LinkedSeqAStarNode(int level, Link link) {
		assert (level <= Short.MAX_VALUE);
		this.level = (short)level;
		this.link = link;
		gscore = Double.NaN;
		hscore = Double.NaN;
	}

	public LinkedSeqAStarNode() {
		this(0, new Link());
	}

	@Override
	public LinkedSeqAStarNode assign(int pos, int rt) {
		return new LinkedSeqAStarNode(level + 1, new Link(link, pos, rt));
	}

	@Override
	public void getAssignments(Assignments assignments) {
		Link link = this.link;
		int i=0;
		while (!link.isRoot()) {
			assignments.assignedPos[i] = link.pos;
			assignments.assignedRTs[i] = link.rt;
			i++;
			link = link.parent;
		}
		assignments.numAssigned = i;
		assignments.sortAssigned();
		assignments.updateUnassigned();
	}

	@Override
	public void getSequence(Sequence seq) {
		Link link = this.link;
		while (!link.isRoot()) {
			SeqSpace.Position pos = seq.seqSpace.positions.get(link.pos);
			SeqSpace.ResType rt = pos.resTypes.get(link.rt);
			seq.set(pos, rt);
			link = link.parent;
		}
	}

	@Override
	public Object getData() {
		return data;
	}

	@Override
	public void setData(Object data) {
		this.data = data;
	}

	public Link getLink() {
		return link;
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
	public String toString() {
		StringBuilder buf = new StringBuilder();
		buf.append("[");
		Link link = this.link;
		while (!link.isRoot()) {
			if (buf.length() > 1) {
				buf.append(", ");
			}
			buf.append(link.pos);
			buf.append(":");
			buf.append(link.rt);
			link = link.parent;
		}
		buf.append("]");
		return buf.toString();
	}
}
