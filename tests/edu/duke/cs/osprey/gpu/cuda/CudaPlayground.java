package edu.duke.cs.osprey.gpu.cuda;

import java.io.IOException;
import java.nio.DoubleBuffer;

import edu.duke.cs.osprey.gpu.cuda.kernels.TestKernel;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;

public class CudaPlayground {
	
	private static final int NumElements = 20;
	private static final int NumRuns = 1000;
	
	public static void main(String[] args)
	throws Exception {
		
		// NOTE: samples and such here:
		// https://github.com/jcuda/jcuda-samples/tree/master/JCudaSamples/src/main/java/jcuda
		
		// info on dynamic parallelism:
		// http://docs.nvidia.com/cuda/cuda-c-programming-guide/#cuda-dynamic-parallelism
		
		// init CUDA
		Context context = new Context(Gpus.get().getGpus().get(0));
		
		hostLoop(context);
		//deviceLoop(context);
	
		context.cleanup();
	}
	
	private static void hostLoop(Context context)
	throws IOException {
		
		TestKernel kernel = new TestKernel(context, NumElements);
		
		// init host buffers
		DoubleBuffer a = kernel.getA();
		DoubleBuffer b = kernel.getB();
		a.clear();
		b.clear();
		for (int i=0; i<NumElements; i++) {
			a.put((double)i);
			b.put((double)i);
		}
		
		Stopwatch stopwatch = new Stopwatch().start();
		
		for (int i=0; i<NumRuns; i++) {
			kernel.uploadAsync();
			kernel.runAsync();
			kernel.downloadSync();
		}
		stopwatch.stop();
		System.out.println("total time: " + TimeFormatter.format(stopwatch.getTimeNs(), 1));
		System.out.println("avg time per op: " + TimeFormatter.format(stopwatch.getTimeNs()/NumRuns, 1));
		
		System.out.println("sum");
		DoubleBuffer out = kernel.getOut();
		out.rewind();
		for (int i=0; i<10; i++) {
			System.out.println(out.get());
		}
		
		kernel.cleanup();
	}
	
	private static void deviceLoop(CUmodule module) {
		
		// get the kernel func
		CUfunction func = new CUfunction();
		JCudaDriver.cuModuleGetFunction(func, module, "loop");
		
		// init host buffers
		float[] a = new float[NumElements];
		float[] b = new float[NumElements];
		for (int i=0; i<NumElements; i++) {
			a[i] = (float)i;
			b[i] = (float)i;
		}
		float[] out = new float[NumElements];
		
		// allocate buffers
		CUdeviceptr pdA = new CUdeviceptr();
		JCudaDriver.cuMemAlloc(pdA, NumElements*Sizeof.FLOAT);
		CUdeviceptr pdB = new CUdeviceptr();
		JCudaDriver.cuMemAlloc(pdB, NumElements*Sizeof.FLOAT);
		CUdeviceptr pdOut = new CUdeviceptr();
		JCudaDriver.cuMemAlloc(pdOut, NumElements*Sizeof.FLOAT);
		
		// prep kernel args
		Pointer pParams = Pointer.to(
			Pointer.to(new int[] {NumElements}),
			Pointer.to(pdA),
			Pointer.to(pdB),
			Pointer.to(pdOut)
		);
		
		Stopwatch stopwatch = new Stopwatch().start();
		
		// upload inputs
		JCudaDriver.cuMemcpyHtoDAsync(pdA, Pointer.to(a), NumElements*Sizeof.FLOAT, null);
		JCudaDriver.cuMemcpyHtoDAsync(pdB, Pointer.to(b), NumElements*Sizeof.FLOAT, null);
		
		// launch the kernel
		JCudaDriver.cuLaunchKernel(
			func,
			1, 1, 1,
			1, 1, 1,
			0,
			null,
			pParams,
			null
		);
		
		// download the output
		JCudaDriver.cuMemcpyDtoH(Pointer.to(out), pdOut, NumElements*Sizeof.FLOAT);
		
		stopwatch.stop();
		System.out.println("total time: " + TimeFormatter.format(stopwatch.getTimeNs(), 1));
		System.out.println("avg time per op: " + TimeFormatter.format(stopwatch.getTimeNs()/NumRuns, 1));
		
		System.out.println("sum");
		for (int i=0; i<10; i++) {
			System.out.println(out[i]);
		}
		
		// cleanup
		JCudaDriver.cuMemFree(pdA);
		JCudaDriver.cuMemFree(pdB);
		JCudaDriver.cuMemFree(pdOut);
	}
}
