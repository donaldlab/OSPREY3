/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests.variationalkstar.continuous;

import Jama.Matrix;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.*;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.function.ToDoubleFunction;

/**
 *
 * @author aditya
 */
public class RKHSFunctionTest {
 
    public static void main (String[] args) { 
	int numDims = 1;
	double[] lBounds = new double[numDims];
	double[] uBounds = new double[numDims];
	for (int i=0; i<uBounds.length; i++) { uBounds[i] = 10; }
	
	double[][] bounds = new double[numDims][2];
	for (int dim=0; dim<bounds.length; dim++) { 
	    bounds[dim][0] = lBounds[dim];
	    bounds[dim][1] = uBounds[dim];
	}
	
	double sigma = 2;
	
	KernelGaussian k = new KernelGaussian(bounds, sigma);
	
	ToDoubleFunction<double[]> f = (point)->(Math.sin(point[0]*2));
	
	RKHSFunction func = new RKHSFunction(k, lBounds, uBounds, f, 8);
	RKHSFunction f2 = 
		new RKHSFunction(
			func.featureMaps,
			func.fitCoeffsIdentity(func.featureMaps, func.referenceFunction));
	
	RKHSFunction f3 = 
		new RKHSFunction(
			func.featureMaps,
			func.fitCoeffsInterpolate(func.featureMaps, func.referenceFunction));
	printFunctionSamples(func, f2, f3, f);

	System.out.println("Integral: "+func.computeIntegral());
    }
    
    static void printFunctionSamples(RKHSFunction f, RKHSFunction f2, RKHSFunction f3, ToDoubleFunction<double[]> r) { 
	double[][] samples = f.gridSample(100, f.domainLB, f.domainUB);
	double[][] mats = new double[samples.length][5];
	for (int i=0; i<samples.length; i++) { 
	    mats[i][4] = f3.eval(samples[i]);
	    mats[i][3] = f2.eval(samples[i]);
	    mats[i][2] = f.eval(samples[i]);
	    mats[i][1] = r.applyAsDouble(samples[i]);
	    mats[i][0] = samples[i][0];
	}
	Matrix m = new Matrix(mats);
	
	PrintWriter writer = null;
	try {
	    writer = new PrintWriter("/home/aditya/Desktop/fvals.csv", "UTF-8");
	    m.print(writer, 5, 5);
	    writer.flush();
	} catch (Exception e) {
	    e.printStackTrace();
	    System.out.println("Printing failed");
	}
    }
}
