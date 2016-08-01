package edu.duke.cs.osprey.partcr.pickers;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

public class WalkingConfPicker implements ConfPicker {

	private int numItersPerConf;
	
	// NOTE: don't use int[] as the key value, primitive arrays doesn't have a good hashCode() method
	private Map<List<Integer>,Integer> pickedConfs;
	
	public WalkingConfPicker() {
		this(1);
	}
	
	public WalkingConfPicker(int numItersPerConf) {
		this.numItersPerConf = numItersPerConf;
		pickedConfs = new HashMap<>(); // yeah, ok to match on instance
	}

	@Override
	public ScoredConf pick(List<ScoredConf> confs) {
		
		for (ScoredConf conf : confs) {
			
			// use Java 8 magic to convert the array to a list
			// sadly, Arrays.toList() won't work here =(
			List<Integer> confList = Arrays.stream(conf.getAssignments()).boxed().collect(Collectors.toList());
			
			// have we picked this conf yet?
			Integer numIters = pickedConfs.get(confList);
			if (numIters == null) {
				
				// nope, let's start
				pickedConfs.put(confList, 1);
				return conf;
				
			} else if (numIters < numItersPerConf) {
				
				// yup, but not enough
				pickedConfs.put(confList, numIters + 1);
				return conf;
			}
			
			// yeah, but picked it too much already, so look to the next conf
		}
		
		// we ran out of confs, start over
		pickedConfs.clear();
		return pick(confs);
	}
}
