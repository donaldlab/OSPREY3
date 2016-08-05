package edu.duke.cs.osprey.tests.variationalkstar.continuous;

import edu.duke.cs.osprey.partitionfunctionbounds.continuous.RKHS_MRF;

public class RKHS_MRFTests {

	public static void main (String[] args) {
		RKHS_MRF mrf = new RKHS_MRF();
		mrf.testEquidistantSampling();
	}

}
