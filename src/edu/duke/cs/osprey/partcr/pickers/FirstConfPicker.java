package edu.duke.cs.osprey.partcr.pickers;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

public class FirstConfPicker implements ConfPicker {

	@Override
	public ScoredConf pick(List<ScoredConf> confs) {
		return confs.get(0);
	}
}
