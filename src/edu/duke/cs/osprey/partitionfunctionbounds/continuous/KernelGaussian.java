package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.Arrays;

import edu.duke.cs.osprey.tools.VectorAlgebra;

public class KernelGaussian extends Kernel {

	private double sigma;
	
	/** 
	 * Gaussian kernel - k(x, y) = exp(-1*||x-y||)
	 * @param bounds
	 */
	public KernelGaussian(double[][] bounds, double sigma) {
		super(bounds);
		this.sigma = sigma;
	}

	// evaluates the kernel at two points in the domain
	public double eval(double[] x, double[] y) {
		if (! super.validInput(x, y)) {
			throw new RuntimeException("Input to Gaussian kernel not valid.");
		}
		return Math.exp(-1 * this.distance(x, y) / this.sigma);
	}
	
	// returns the Euclidean distance between x and y
	double distance(double[] x, double[] y) {
		double dist = 0.0;
		for (int i=0; i<x.length; i++) {
			dist += Math.pow((x[i] - y[i]), 2);
		}
		return Math.sqrt(dist);
	}

}
