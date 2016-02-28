package edu.duke.cs.osprey.kstar.implementation;

import java.util.HashMap;

import edu.duke.cs.osprey.confspace.AllowedSeqs;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.KSAbstract;

public class KSImplSubLinear extends KSAbstract {

	public KSImplSubLinear(ConfigFileParser cfp) {
		super(cfp);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void init(HashMap<Integer, AllowedSeqs> strand2AllowedSeqs) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String getKSMethod() {
		return "sublinear";
	}

}
