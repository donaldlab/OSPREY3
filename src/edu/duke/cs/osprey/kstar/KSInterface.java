package edu.duke.cs.osprey.kstar;

import java.util.HashMap;

import edu.duke.cs.osprey.confspace.AllowedSeqs;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public interface KSInterface {

	public void init( HashMap<Integer, AllowedSeqs> strand2AllowedSeqs );
	
	public void run();
	
	public String getKSMethod();
	
}
