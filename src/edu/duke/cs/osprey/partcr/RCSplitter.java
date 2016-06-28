package edu.duke.cs.osprey.partcr;

import java.util.List;

import edu.duke.cs.osprey.confspace.RC;

public interface RCSplitter {

	boolean willSplit(int pos, RC rc);
	List<RC> split(int pos, RC rc);
}
