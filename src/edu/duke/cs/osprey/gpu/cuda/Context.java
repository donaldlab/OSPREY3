package edu.duke.cs.osprey.gpu.cuda;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

import edu.duke.cs.osprey.tools.ResourceExtractor;
import jcuda.Pointer;
import jcuda.driver.CUcontext;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;

public class Context {
	
	private Gpu gpu;
	private CUcontext context;
	
	private Map<String,CUmodule> kernels;
	
	public Context(Gpu gpu) {
		
		this.gpu = gpu;
		
		// create the cuda context
		context = new CUcontext();
		JCudaDriver.cuCtxCreate(context, 0, gpu.getDevice());
		
		kernels = new HashMap<>();
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	public synchronized CUmodule getKernel(String name)
	throws IOException {
		
		// check the cache first
		CUmodule kernel = kernels.get(name);
		if (kernel == null) {
				
			// cache miss, load the kernel
		
			// do we have a kernel binary?
			String resourcePath = String.format("kernelBinaries/%s.bin", name);
			URL url = getClass().getResource(resourcePath);
			if (url == null) {
				throw new IOException("precompiled kernel binary not found at: " + resourcePath);
			}
			
			// yup, load it
			File kernelFile = ResourceExtractor.extract(url);
			kernel = new CUmodule();
			JCudaDriver.cuModuleLoad(kernel, kernelFile.getAbsolutePath());
			
			// update the cache
			kernels.put(name, kernel);
		}
		
		return kernel;
	}
	
	public CUdeviceptr malloc(long numBytes) {
		CUdeviceptr pdBuf = new CUdeviceptr();
		JCudaDriver.cuMemAlloc(pdBuf, numBytes);
		return pdBuf;
	}
	
	public void free(CUdeviceptr pdBuf) {
		JCudaDriver.cuMemFree(pdBuf);
	}
	
	public void uploadSync(CUdeviceptr pdBuf, Pointer phBuf, long numBytes) {
		JCudaDriver.cuMemcpyHtoD(pdBuf, phBuf, numBytes);
	}
	
	public void uploadAsync(CUdeviceptr pdBuf, Pointer phBuf, long numBytes, GpuStream stream) {
		JCudaDriver.cuMemcpyHtoDAsync(pdBuf, phBuf, numBytes, stream.getStream());
	}
	
	public void downloadSync(Pointer phBuf, CUdeviceptr pdBuf, long numBytes) {
		JCudaDriver.cuMemcpyDtoH(phBuf, pdBuf, numBytes);
	}
	
	public void launchKernel(CUfunction func, int gridBlocks, int blockThreads, int sharedMemBytes, Pointer pArgs, GpuStream stream) {
		JCudaDriver.cuLaunchKernel(
			func,
			gridBlocks, 1, 1,
			blockThreads, 1, 1,
			sharedMemBytes,
			stream.getStream(),
			pArgs,
			null
		);
	}
	
	public void waitForGpu() {
		JCudaDriver.cuCtxSynchronize();
	}
	
	public void attachCurrentThread() {
		JCudaDriver.cuCtxPushCurrent(context);
	}
	
	public void detatchCurrentThread() {
		JCudaDriver.cuCtxPopCurrent(context);
	}
	
	public synchronized void cleanup() {
		
		for (CUmodule kernel : kernels.values()) {
			JCudaDriver.cuModuleUnload(kernel);
		}
		kernels.clear();
		
		JCudaDriver.cuCtxDestroy(context);
	}
}
