package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public interface KSInterface {

	public void init( HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs );
	
	public void createEmats( ArrayList<Boolean> contSCFlexVals );
	
	public void run();
	
	public String getKSMethod();
	
}
