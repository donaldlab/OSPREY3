package edu.duke.cs.osprey.minimization;

public class MinimizerFactory {

	private static String implementation = "ccd";
	
	public static void setImplementation( String impl ) {
		implementation = impl;
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
