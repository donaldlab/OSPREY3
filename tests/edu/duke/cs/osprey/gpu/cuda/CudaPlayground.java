package edu.duke.cs.osprey.gpu.cuda;

import java.io.File;

import edu.duke.cs.osprey.tools.ResourceExtractor;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUcontext;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;

public class CudaPlayground {
	
	private static final int NumElements = 10000;
	private static final int NumRuns = 10;
	
	private static final int BlockThreads = 256;
	private static final int GridBlocks = 40;//(NumElements + BlockThreads - 1)/BlockThreads;
	
	public static void main(String[] args)
	throws Exception {
		
		// NOTE: samples and such here:
		// https://github.com/jcuda/jcuda-samples/tree/master/JCudaSamples/src/main/java/jcuda
		
		// info on dynamic parallelism:
		// http://docs.nvidia.com/cuda/cuda-c-programming-guide/#cuda-dynamic-parallelism
		
		// init CUDA
		Gpu gpu = Gpus.get().getGpus().get(0);
		CUcontext context = gpu.makeContextForThisThread();
		
		// load the kernel
		File kernelFile = ResourceExtractor.extract("kernelBinaries/test.cubin", Gpu.class);
		CUmodule module = new CUmodule();
		JCudaDriver.cuModuleLoad(module, kernelFile.getAbsolutePath());
		
		hostLoop(module);
		deviceLoop(module);
		
		JCudaDriver.cuModuleUnload(module);
		JCudaDriver.cuCtxDestroy(context);
	}
	
	private static void hostLoop(CUmodule module) {
		
		// get the kernel func
		CUfunction func = new CUfunction();
		JCudaDriver.cuModuleGetFunction(func, module, "add");
		
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
		
		for (int i=0; i<NumRuns; i++) {
		
			// upload inputs
			JCudaDriver.cuMemcpyHtoDAsync(pdA, Pointer.to(a), NumElements*Sizeof.FLOAT, null);
			JCudaDriver.cuMemcpyHtoDAsync(pdB, Pointer.to(b), NumElements*Sizeof.FLOAT, null);
			
			// launch the kernel
			JCudaDriver.cuLaunchKernel(
				func,
				GridBlocks, 1, 1,
				BlockThreads, 1, 1,
				0,
				null,
				pParams,
				null
			);
			
			// download the output
			JCudaDriver.cuMemcpyDtoH(Pointer.to(out), pdOut, NumElements*Sizeof.FLOAT);
		}
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
