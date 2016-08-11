package edu.duke.cs.osprey.partcr.pickers;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

public interface ConfPicker {
	
	ScoredConf pick(List<ScoredConf> confs);
}
