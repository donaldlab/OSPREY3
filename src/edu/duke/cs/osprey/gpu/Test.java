package edu.duke.cs.osprey.gpu;

import java.nio.DoubleBuffer;
import java.util.Random;

import org.apache.commons.io.IOUtils;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLMemory;
import com.jogamp.opencl.CLPlatform;
import com.jogamp.opencl.CLProgram;

public class Test {
	public static void main(String[] args)
	throws Exception {
		
		// init random doubles
		System.out.println("creating data...");
		final int n = 1024 * 1024;
		double[] a = new double[n];
		double[] b = new double[n];
		double[] out = new double[n];
		Random rand = new Random(12345);
		for (int i=0; i<n; i++) {
			a[i] = rand.nextDouble();
			b[i] = rand.nextDouble();
		}
		
		final int NumRuns = 1;
		
		// benchmark vector addition on cpu
		System.out.println("benchmarking CPU...");
		long timeNs = System.nanoTime();
		for (int i=0; i<NumRuns; i++) {
			addCpu(a, b, out);
		}
		timeNs = System.nanoTime() - timeNs;
		System.out.println(String.format("time: %d ms", timeNs/1000000));
		
		// check opencl platforms
		System.out.println("platforms:");
		for (CLPlatform platform : CLPlatform.listCLPlatforms()) {
			System.out.println("\t" + platform);
			for (String extension : platform.getExtensions()) {
				System.out.println("\t\t" + extension);
			}
		}
		
		// init opencl
		CLContext context = CLContext.create();
		try {
			
			// init the gpu
			System.out.println("platform: " + context.getPlatform());
			CLDevice gpu = context.getMaxFlopsDevice(CLDevice.Type.GPU);
			System.out.println("gpu: " + gpu);
			if (gpu == null) {
				throw new Error("no GPU!");
			}
			CLCommandQueue commandQueue = gpu.createCommandQueue();
			
			// create transfer buffers
			int groupSize = gpu.getMaxWorkGroupSize();
			int workSize = n;
			workSize = roundUpWorkSize(workSize, groupSize);
			CLBuffer<DoubleBuffer> bufA = context.createDoubleBuffer(workSize, CLMemory.Mem.READ_ONLY);
			CLBuffer<DoubleBuffer> bufB = context.createDoubleBuffer(workSize, CLMemory.Mem.READ_ONLY);
			CLBuffer<DoubleBuffer> bufOut = context.createDoubleBuffer(workSize, CLMemory.Mem.WRITE_ONLY);
			
			// copy data to buffers
			for (int i=0; i<n; i++) {
				bufA.getBuffer().put((float)a[i]);
				bufB.getBuffer().put((float)b[i]);
			}
			bufA.getBuffer().rewind();
			bufB.getBuffer().rewind();
			
			// load the kernel
			CLProgram program = context.createProgram(IOUtils.toString(Test.class.getResourceAsStream("kernels/add.cl"))).build();
			CLKernel kernel = program.createCLKernel("addDouble");
			kernel.putArgs(bufA, bufB, bufOut);
			
			// upload args to gpu (nonblocking)
			commandQueue.putWriteBuffer(bufA, false);
			commandQueue.putWriteBuffer(bufB, false);
			
			System.out.println("benchmarking GPU...");
			timeNs = System.nanoTime();
			for (int i=0; i<NumRuns; i++) {
				
				// run kernel, download results (blocking)
				commandQueue.put1DRangeKernel(kernel, 0, workSize, groupSize);
				commandQueue.putReadBuffer(bufOut, true);
			}
			timeNs = System.nanoTime() - timeNs;
			System.out.println(String.format("time: %d ms", timeNs/1000000));
			
			// TEMP: check accuracy
			for (int i=0; i<10; i++) {
				System.out.println(String.format("a: %f, b: %f, cpu: %f, gpu: %f",
					a[i], b[i], out[i], bufOut.getBuffer().get(i)
				));
			}
			
		} finally {
			context.release();
		}
		
		System.out.println("done!");
	}
	
	private static void addCpu(double[] a, double[] b, double[] out) {
		int n = out.length;
		for (int i=0; i<n; i++) {
			out[i] = a[i] + b[i];
		}
	}
	
	private static int roundUpWorkSize(int workSize, int groupSize) {
		int r = groupSize % workSize;
		if (r == 0) {
			return groupSize;
		} else {
			return groupSize + workSize - r;
		}
	}
}
