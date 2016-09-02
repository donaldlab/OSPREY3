package edu.duke.cs.osprey.tests.variationalkstar.continuous;

import java.util.function.ToDoubleFunction;

import edu.duke.cs.osprey.partitionfunctionbounds.continuous.KernelGaussian;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.RKHSFunction;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.RKHS_MRF;

public class RKHS_MRFTests {

	public static void main (String[] args) {
		TestFunction tf1 = new TestFunction(1.00);
		TestFunction tf2 = new TestFunction(2.00);
		
		double[] upperBounds = new double[]{0};
		double[] lowerBounds = new double[]{180};
		double[][] bounds = new double[][]{upperBounds, lowerBounds};
		KernelGaussian k = new KernelGaussian(bounds, /*sigma*/ 15);

		RKHSFunction rf1 = new RKHSFunction(k, lowerBounds, upperBounds, tf1);
		RKHSFunction rf2 = new RKHSFunction(k, lowerBounds, upperBounds, tf2);
		
		MixFunction mf = new MixFunction(tf1, tf2);
		
		RKHSFunction rfm = new RKHSFunction(k, lowerBounds, upperBounds, mf);
	}
	
	private static class MixFunction implements ToDoubleFunction<double[]> {
		private TestFunction tf1;
		private TestFunction tf2;
		
		public MixFunction(TestFunction f1, TestFunction f2) {
			this.tf1 = f1;
			this.tf2 = f2;
		}
		
		public double applyAsDouble(double[] value) { 
			return this.tf1.applyAsDouble(value) * this.tf2.applyAsDouble(value);
		}
	}
	
	private static class TestFunction implements ToDoubleFunction<double[]> {
		
		private double sigma;
		
		public TestFunction(double s) {
			sigma = s;
		}
		@Override
		public double applyAsDouble(double[] value) {
			// TODO Auto-generated method stub
			double result = 1.00;
			for (double v : value) { 
				result += 10;
			}
			return result;
		}
		
	}

}
