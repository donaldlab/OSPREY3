package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

public abstract class Kernel {
	
	/**
	 * A kernel on a domain D is a symmetric positive-definite function 
	 * 		k: D x D -> R
	 * Here the domain is assumed to be a subset of some n-dimensional 
	 * Euclidean space, i.e. [0, 1]^n.
	 * 
	 * Specific kernels should be in classes named along the lines of 
	 * KernelGaussian or KernelPoisson
	 */
	
	// bounds[i][0] --> lower bound on the ith dimension
	// bounds[i][1] --> upper bound on the ith dimension
	double[][] bounds; 

	public Kernel(double[][] bounds) {
		this.bounds = bounds;
	}
	
	public abstract double eval(double[] x, double[] y);
	
	public boolean validInput(double[] x, double[] y) {
		boolean isValid = true;
		if (x.length != y.length) { isValid = false; }
		for (int i=1; i<x.length; i++) {
			boolean xLB = (x[i] > bounds[i][0]); 
			boolean xUB = (x[i] < bounds[i][1]);
			boolean yLB = (y[i] > bounds[i][0]);
			boolean yUB = (y[i] < bounds[i][1]);
			isValid = isValid && xLB && xUB && yLB && yUB;
		}
		return isValid;
	}
}
