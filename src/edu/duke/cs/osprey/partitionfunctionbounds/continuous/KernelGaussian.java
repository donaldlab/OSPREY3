package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.Arrays;

import edu.duke.cs.osprey.tools.VectorAlgebra;

/**
 * Implements a Gaussian kernel.
 * @author aditya
 */
public class KernelGaussian extends Kernel {
    
    // Gaussian kernel has a variance (sigma) and associated bounds 
    private double sigma;
    private double[][] bounds;
    
    public KernelGaussian(double[][] domainBounds, double sigma) {
        super(domainBounds);
        this.bounds = domainBounds;
        this.sigma = sigma;
    }
    
    /**
     * Evaluates the kernel at a specific pair of points 
     * @param x
     * @param y
     * @return 
     */
    public double eval(double[] x, double[] y) {
        if (! super.validInput(x, y)) {
            throw new RuntimeException("Input to Gaussian kernel not valid: " +
                    "input was:\n\t"+
		    Arrays.toString(x)+  
		    "\n\t"+Arrays.toString(y)+
		    "\nbut bounds were\n"+
		    printBounds());
        }
        return Math.exp(-1 * (Math.pow(this.distance(x, y),2)/ (2 * Math.pow(this.sigma,2))));
    }
    
    /**
     * Returns the Euclidean distance between two points 
     * @param x
     * @param y
     * @return 
     */
    double distance(double[] x, double[] y) {
        double dist = 0.0;
        for (int i=0; i<x.length; i++) {
            dist += Math.pow((x[i] - y[i]), 2);
        }
        return Math.sqrt(dist);
    }
    
    /**
     * Returns a string representation of the bounds of the kernel domain 
     * @return 
     */
    String printBounds() {
        StringBuilder s = new StringBuilder();
	for (int i=0; i<bounds.length; i++) { 
	    double[] bound = bounds[i];
	    s.append("\tdim " + i + ": "+Arrays.toString(bound) + "\n");
	}
        return s.toString();
    }
    
}
