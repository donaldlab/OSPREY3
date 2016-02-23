package edu.duke.cs.osprey.kstar;

import java.util.HashMap;

import edu.duke.cs.osprey.confspace.AllowedSeqs;

public interface KSInterface {

	public void init( HashMap<Integer, AllowedSeqs> strand2AllowedSeqs );
	
	public void run();
	
	public String getKSMethod();
	
}
