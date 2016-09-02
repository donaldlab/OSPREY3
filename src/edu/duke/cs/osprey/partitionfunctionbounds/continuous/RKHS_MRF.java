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
		int numSamples = 2000;
		DoubleFunction<Double> f = (double x) -> (0.1);		
		
		double[] domainLB = new double[]{0};
		double[] domainUB = new double[]{10};
		double[][] domainBounds = new double[][]{domainLB, domainUB};
		
		k = new KernelGaussian(domainBounds, 1);
		
		double[] equiSamples = sampleDomain(numSamples, domainLB[0], domainUB[0]);
		FeatureMap[] featureMaps = getFeatureMaps(equiSamples, k);
		double[] coeffs = fitCoeffs(featureMaps, f, k);
		
		double[] testSamples = sampleDomain(numSamples*10, domainLB[0], domainUB[0]);
		double rmsd=0.0;
		for (double d : testSamples) {
			double delta = f.apply(d) - evalLC(featureMaps, coeffs, new double[]{d});
			rmsd += Math.pow(delta, 2);
		}
		System.out.println("RMSD: "+Math.sqrt(rmsd));
	}
	
	public double[] sampleDomain(int numSamples, double lb, double ub) {
		double[] samples = new double[numSamples];
		int ind = 0;
		double increment = ((ub-lb)/numSamples) * 0.95;
		double i = lb;
		while (ind < numSamples) {
			samples[ind] = i;
			i += increment;
			ind++;
		}
		return samples;
	}

	public FeatureMap[] getFeatureMaps(double[] samples, Kernel k) {
		FeatureMap[] fmaps = new FeatureMap[samples.length];
		for (int i=0; i<fmaps.length; i++) {
			fmaps[i] = new FeatureMap(k, new double[]{samples[i]});
		}
		return fmaps;
	}

	public double[] fitCoeffs (FeatureMap[] fMaps, DoubleFunction<Double> func, Kernel k) {
		Matrix gMat = new Matrix(fMaps.length, fMaps.length);
		for (int i=0; i<fMaps.length; i++) {
			for (int j=0; j<fMaps.length; j++) {
				gMat.set(i, j, k.eval(fMaps[i].getLoc(), fMaps[j].getLoc()));
			}
		}
		
		Matrix fMat = new Matrix(fMaps.length, 1);
		for (int i=0; i<fMaps.length; i++) {
			fMat.set(i, 0, func.apply(fMaps[i].getLoc()[0]));
		}
		
		Matrix cMat = gMat.solve(fMat);
		return cMat.transpose().getArray()[0];
	}

	public double evalLC (FeatureMap[] fMaps, double[] coeffs, double[] x) {
		double result = 0.0;
		for (int i=0; i<fMaps.length; i++) {
			result += fMaps[i].eval(x) * coeffs[i];
		}
		return result;
	}
}
