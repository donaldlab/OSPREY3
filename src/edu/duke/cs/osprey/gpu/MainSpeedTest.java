package edu.duke.cs.osprey.gpu;

import java.util.Random;

import edu.duke.cs.osprey.gpu.kernels.TestAddKernel;
import edu.duke.cs.osprey.gpu.kernels.TestFancyKernel;

@SuppressWarnings("unused")
public class MainSpeedTest {
	
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
		}
		
		final int NumRuns = 10;
		
		// benchmark vector addition on cpu
		System.out.println("benchmarking CPU...");
		long cpuTimeNs = System.nanoTime();
		for (int i=0; i<NumRuns; i++) {
			//cpuAdd(a, b, out);
			cpuFancy(a, b, out);
		}
		cpuTimeNs = System.nanoTime() - cpuTimeNs;
		System.out.println(String.format("time: %d ms", cpuTimeNs/1000000));
		
		// init the kernel
		//TestAddKernel.Bound kernel = new TestAddKernel().bind();
		TestFancyKernel.Bound kernel = new TestFancyKernel().bind();
			
		// copy data to buffers
		kernel.setArgs(n);
		for (int i=0; i<n; i++) {
			kernel.getA().put(a[i]);
			kernel.getB().put(b[i]);
		}
		kernel.getA().rewind();
		kernel.getB().rewind();
			
		// upload args to gpu
		kernel.uploadSync();
			
		System.out.println("benchmarking GPU...");
		long gpuTimeNs = System.nanoTime();
		for (int i=0; i<NumRuns; i++) {
			kernel.runAsync();
		}
		kernel.waitForGpu();
		gpuTimeNs = System.nanoTime() - gpuTimeNs;
		System.out.println(String.format("time: %d ms => %.2fx speedup", gpuTimeNs/1000000, (float)cpuTimeNs/gpuTimeNs));
			
		// TEMP: check accuracy
		kernel.downloadSync();
		for (int i=0; i<10; i++) {
			double gpuVal = kernel.getOut().get(i);
			System.out.println(String.format("a: %f, b: %f, cpu: %f, gpu: %f",
				a[i], b[i], out[i], gpuVal
			));
			double err = Math.abs(out[i] - gpuVal);
			assert (err == 0) : "err: " + err;
		}
		
		System.out.println("done!");
	}
	
	private static void cpuAdd(double[] a, double[] b, double[] out) {
		int n = out.length;
		for (int i=0; i<n; i++) {
			out[i] = a[i] + b[i];
		}
	}
	
	private static void cpuFancy(double[] a, double[] b, double[] out) {
		int n = out.length;
		for (int i=0; i<n; i++) {
			out[i] = Math.sqrt(a[i]*a[i] + b[i]*b[i]);
		}
	}
}
