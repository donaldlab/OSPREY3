package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.Random;
import java.util.function.DoubleFunction;

import Jama.*;

public class RKHS_MRF {

	Kernel k;
	
	public RKHS_MRF() {
		// TODO Auto-generated constructor stub

	}
	
	public void testEquidistantSampling() {
		// set up domain as [0, 10]^numDims
		int numDims = 1;
		double[][] bounds = new double[numDims][2];
		for (int i=0; i<bounds.length; i++) {
			bounds[i][0] = 0;
			bounds[i][1] = 10;
		}

		// set up kernel
		this.k = new KernelGaussian(bounds);
		// f = func(x, 1); g = func(x, 2);

		// set up linear combination of equidistant samples
		FeatureMap[] fMaps = getFeatureMaps(9);
		
		// compute optimal coefficients
		// 		first, get Gram matrix
		Matrix gMat = new Matrix(9, 9);
		for (int i=0; i<9; i++) {
			for (int j=0; j<9; j++) {
				gMat.set(i, j, k.eval(toArr(i+1), toArr(j+1))); 	
			}
		}
		// 		then, get vector of function values
		Matrix fMat = new Matrix(9, 1);
		for (int i=0; i<9; i++) {
			fMat.set(i, 0, func(toArr(i+1), 1));
		}
		// 		then, solve for x in gMat*X = fMat
		Matrix cMat = gMat.solve(fMat);
		double[] coeffs = new double[9];
		for (int i=0; i<coeffs.length; i++) { coeffs[i] = cMat.get(i, 0); }
		
		DoubleFunction<Double> fitted = (double x) -> evalLC(fMaps, coeffs, toArr(x));

		double rmsd = 0.0;
		for (double i=0.0; i<10; i+=0.5) {
			double funcVal = func(toArr(i), 1);
			double fitVal = fitted.apply(i);
			double error = Math.sqrt((funcVal-fitVal)*(funcVal-fitVal));
			System.out.println("At "+i+" error was "+(error/funcVal * 100)+"%");
			System.out.println("    Value: "+funcVal);
			System.out.println("    Prediction: "+fitVal);
		}
		System.out.println("RMSD: "+Math.sqrt(rmsd));
	}
	
	public double evalLC(FeatureMap[] fMaps, double[] coeffs, double[] point) {
		double val = 0;
		for (int i=0; i<fMaps.length; i++) {
			val += fMaps[i].eval(point)*coeffs[i];
		}
		return val;
	}
	
	public double[] toArr(double x) {
		double[] a = new double[1];
		a[0] = x;
		return a;
	}
	
	public FeatureMap[] getFeatureMaps(int n) {
		FeatureMap[] fmaps = new FeatureMap[n];
		for (int i=0; i<fmaps.length; i++) {
			double[] curPoint = new double[1];
			curPoint[0] = i;
			fmaps[i] = new FeatureMap(this.k, curPoint);
		}
		return fmaps;
	}
	
	public double func(double[] x, int seed) {
		Random random = new Random(seed);
		int min = 0;
		int max = 100;
		
		// set up array of feature maps
		FeatureMap[] fmaps = getFeatureMaps(9);
		
		// set up array of linear coefficients
		double[] coeffs = new double[fmaps.length];
		for (int i=0; i<coeffs.length; i++) {
			coeffs[i] = random.nextDouble() * (max-min) + min;
		}
		
		// take dot product
		double val = 0.0;
		for (int i=0; i<fmaps.length; i++) {
			val += fmaps[i].eval(x)*coeffs[i];
		}
		
//		return val;
		return 30;
	}

}
