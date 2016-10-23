package edu.duke.cs.osprey.gpu.cuda;

import java.io.File;

import edu.duke.cs.osprey.tools.ResourceExtractor;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUcontext;
import jcuda.driver.CUdevice;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;

public class CudaPlayground {
	
	public static void main(String[] args)
	throws Exception {
		
		// NOTE: samples and such here:
		// https://github.com/jcuda/jcuda-samples/tree/master/JCudaSamples/src/main/java/jcuda
		
		// info on dynamic parallelism:
		// http://docs.nvidia.com/cuda/cuda-c-programming-guide/#cuda-dynamic-parallelism
		
		JCudaDriver.setExceptionsEnabled(true);
		
		// init CUDA
		JCudaDriver.cuInit(0);
		CUdevice device = new CUdevice();
		JCudaDriver.cuDeviceGet(device, 0);
		CUcontext context = new CUcontext();
		JCudaDriver.cuCtxCreate(context, 0, device);
		
		// load the kernel
		File kernelFile = ResourceExtractor.extract("kernelSource/test.ptx", Kernel.class);
		CUmodule module = new CUmodule();
		JCudaDriver.cuModuleLoad(module, kernelFile.getAbsolutePath());
		
		// get the kernel func
		CUfunction func = new CUfunction();
		JCudaDriver.cuModuleGetFunction(func, module, "add");
		
		// init host buffers
		final int NumElements = 2000;
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
		int blockSizeX = 256;
		int gridSizeX = (NumElements + blockSizeX - 1)/blockSizeX;
		
		final int NumRuns = 1000;
		//Profiler p = new Profiler();
		Stopwatch stopwatch = new Stopwatch().start();
		
		for (int i=0; i<NumRuns; i++) {
		
			//p.start("upload");
		
			// upload inputs
			//JCudaDriver.cuMemcpyHtoDAsync(pdA, Pointer.to(a), NumElements*Sizeof.FLOAT, stream);
			JCudaDriver.cuMemcpyHtoD(pdA, Pointer.to(a), NumElements*Sizeof.FLOAT);
			JCudaDriver.cuMemcpyHtoD(pdB, Pointer.to(b), NumElements*Sizeof.FLOAT);
			
			//p.start("kernel");
			
			// launch the kernel
			JCudaDriver.cuLaunchKernel(
				func,
				gridSizeX, 1, 1,
				blockSizeX, 1, 1,
				0,
				null,
				pParams,
				null
			);
			JCudaDriver.cuCtxSynchronize();
			
			//p.start("download");
			
			// download the output
			JCudaDriver.cuMemcpyDtoH(Pointer.to(out), pdOut, NumElements*Sizeof.FLOAT);
			
			//p.stop();
		}
		//System.out.println(p.makeReport(TimeUnit.MICROSECONDS));
		stopwatch.stop();
		System.out.println("avg time per op: " + TimeFormatter.format(stopwatch.getTimeNs()/NumRuns, 1));
		
		System.out.println("sum");
		for (int i=0; i<10; i++) {
			System.out.println(out[i]);
		}
		
		// cleanup
		JCudaDriver.cuMemFree(pdA);
		JCudaDriver.cuMemFree(pdB);
		JCudaDriver.cuMemFree(pdOut);
		JCudaDriver.cuModuleUnload(module);
	}
}
