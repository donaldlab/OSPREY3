package edu.duke.cs.osprey.gpu;

import java.util.Random;

import edu.duke.cs.osprey.gpu.kernels.TestFancyKernel;

public class MainAccuracyTest {
	
	public static void main(String[] args)
	throws Exception {
		
		// init random doubles
		System.out.println("creating data...");
		final int n = 1024 * 1024 * 16;
		double[] a = new double[n];
		double[] b = new double[n];
		double[] out = new double[n];
		Random rand = new Random(12345);
		for (int i=0; i<n; i++) {
			a[i] = rand.nextDouble();
			b[i] = rand.nextDouble();
			out[i] = Math.sqrt(a[i]*a[i] + b[i]*b[i]);
		}
		
		final int NumRuns = 1000;
		
		TestFancyKernel.Bound kernel = new TestFancyKernel().bind();
		
		// copy data to buffers
		kernel.setWorkSize(n);
		for (int i=0; i<n; i++) {
			kernel.getA().put(a[i]);
			kernel.getB().put(b[i]);
		}
		kernel.getA().rewind();
		kernel.getB().rewind();
		
		// upload args to gpu
		kernel.uploadSync();
		
		System.out.println("checking gpu...");
		for (int i=0; i<NumRuns; i++) {
			kernel.runSync();
			kernel.downloadSync();
			for (int j=0; j<n; j++) {
				double gpuVal = kernel.getOut().get(j);
				double err = Math.abs(out[j] - gpuVal);
				assert (err <= 1e-15) : "err: " + err;
			}
		}
		
		System.out.println("done, no errors detected");
	}
}
