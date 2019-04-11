package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class DofInfo {

	public final RCTuple tuple;

	public final List<Integer> counts = new ArrayList<>();
	public final List<Integer> offsets = new ArrayList<>();
	public final List<Integer> ids = new ArrayList<>();
	public final List<Strand> strands = new ArrayList<>();
	public final List<SimpleConfSpace.Position> positions = new ArrayList<>();
	public final List<SimpleConfSpace.ResidueConf> resConfs = new ArrayList<>();

	public final Map<String,Integer> blockIndicesByResNum = new HashMap<>();

	private int numDofs = 0;

	public DofInfo(RCTuple tuple) {
		this.tuple = tuple;
	}

	public int addStrand(Strand strand, int numDofs) {
		int blockIndex = counts.size();
		counts.add(numDofs);
		offsets.add(numDofs());
		ids.add(-1); // TODO: what to use here?
		strands.add(strand);
		positions.add(null);
		resConfs.add(null);
		this.numDofs += numDofs;
		return blockIndex;
	}

	public int addPos(SimpleConfSpace.Position pos, SimpleConfSpace.ResidueConf rc, int numDofs) {
		int blockIndex = counts.size();
		counts.add(numDofs);
		offsets.add(numDofs());
		ids.add(pos.index);
		strands.add(null);
		positions.add(pos);
		resConfs.add(rc);
		blockIndicesByResNum.put(pos.resNum, blockIndex);
		this.numDofs += numDofs;
		return blockIndex;
	}

	public int size() {
		return counts.size();
	}

	public int numDofs() {
		return numDofs;
	}

	public Integer getBlockIndex(String resNum) {
		return blockIndicesByResNum.get(resNum);
	}

	@Override
	public String toString() {
		StringBuilder buf = new StringBuilder();
		for (int i=0; i<size(); i++) {
			buf.append(String.format("block %s:\n", i));
			buf.append(String.format("\toffset:      %d\n", offsets.get(i)));
			buf.append(String.format("\tcount:       %d\n", counts.get(i)));
			buf.append(String.format("\tstrand:      %s\n", strands.get(i)));
			buf.append(String.format("\tpos:         %d %s\n", positions.get(i).index, positions.get(i).resNum));
			buf.append(String.format("\trc:          %d %s %s\n", resConfs.get(i).index, resConfs.get(i).template.name, resConfs.get(i).getRotamerCode()));
		}
		return buf.toString();
	}
}
