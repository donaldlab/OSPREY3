package edu.duke.cs.osprey.partcr.pickers;

import java.util.List;

public class FirstConfPicker implements ConfPicker {

	@Override
	public int[] pick(List<int[]> confs) {
		return confs.get(0);
	}
}
