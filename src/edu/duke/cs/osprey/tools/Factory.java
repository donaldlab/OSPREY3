package edu.duke.cs.osprey.tools;

public interface Factory<VT,CT> {
	
	VT make(CT context);
}
