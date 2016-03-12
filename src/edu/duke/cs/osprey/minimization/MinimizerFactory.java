package edu.duke.cs.osprey.minimization;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class MinimizerFactory {

	private static String implementation = "ccd";

	public static void setImpl( String impl ) {

		switch( impl ) {
		
		case "bfgs":
		case "hbfgsccd":
		case "ccd":
			implementation = impl;
			return;

		default:
			throw new RuntimeException("ERROR: invalid implementation.");
		}


	}


	public static String getImpl() {
		return implementation;
	}


	public static Minimizer getMinimizer( ObjectiveFunction ofn, boolean useCorners ) {

		switch( implementation ) {

		case "bfgs":
			return new BFGSMinimizer(ofn, useCorners);

		case "hbfgsccd":
			return new HBFGSCCDMinimizer(ofn, useCorners);

		case "ccd":
		default:
			return new CCDMinimizer(ofn, useCorners);

		}

	}

}
