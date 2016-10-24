package edu.duke.cs.osprey.gpu.cuda;

import java.io.File;
import java.io.IOException;
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
	private Thread thread;
	
	private Map<String,CUmodule> kernels;
	
	public Context(Gpu gpu) {
		
		this.gpu = gpu;
		
		// create the cuda context
		context = new CUcontext();
		JCudaDriver.cuCtxCreate(context, 0, gpu.getDevice());
		
		// keep track of which thread created this context
		// only access by this thread is allowed in the future
		thread = Thread.currentThread();
		
		kernels = new HashMap<>();
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	public CUmodule getKernel(String kernelFilename)
	throws IOException {
		checkThread();
		
		// check the cache first
		CUmodule kernel = kernels.get(kernelFilename);
		if (kernel == null) {
		
			// cache miss, load the kernel
			File kernelFile = ResourceExtractor.extract("kernelBinaries/" + kernelFilename, getClass());
			kernel = new CUmodule();
			JCudaDriver.cuModuleLoad(kernel, kernelFile.getAbsolutePath());
			
			// update the cache
			kernels.put(kernelFilename, kernel);
		}
		
		return kernel;
	}
	
	public CUdeviceptr malloc(long numBytes) {
		checkThread();
		CUdeviceptr pdBuf = new CUdeviceptr();
		JCudaDriver.cuMemAlloc(pdBuf, numBytes);
		return pdBuf;
	}
	
	public void free(CUdeviceptr pdBuf) {
		checkThread();
		JCudaDriver.cuMemFree(pdBuf);
	}
	
	public void uploadSync(CUdeviceptr pdBuf, Pointer phBuf, long numBytes) {
		checkThread();
		JCudaDriver.cuMemcpyHtoD(pdBuf, phBuf, numBytes);
	}
	
	public void uploadAsync(CUdeviceptr pdBuf, Pointer phBuf, long numBytes) {
		checkThread();
		JCudaDriver.cuMemcpyHtoDAsync(pdBuf, phBuf, numBytes, null);
	}
	
	public void downloadSync(Pointer phBuf, CUdeviceptr pdBuf, long numBytes) {
		checkThread();
		JCudaDriver.cuMemcpyDtoH(phBuf, pdBuf, numBytes);
	}
	
	public void launchKernel(CUfunction func, int gridBlocks, int blockThreads, int sharedMemBytes, Pointer pArgs) {
		checkThread();
		JCudaDriver.cuLaunchKernel(
			func,
			gridBlocks, 1, 1,
			blockThreads, 1, 1,
			sharedMemBytes,
			null,
			pArgs,
			null
		);
	}
	
	public void waitForGpu() {
		checkThread();
		JCudaDriver.cuCtxSynchronize();
	}
	
	public void cleanup() {
		checkThread();
		
		for (CUmodule kernel : kernels.values()) {
			JCudaDriver.cuModuleUnload(kernel);
		}
		kernels.clear();
		
		JCudaDriver.cuCtxDestroy(context);
	}
	
	private void checkThread() {
		Thread currentThread = Thread.currentThread();
		if (!currentThread.equals(thread)) {
			throw new Error(String.format("CUDA context created in thread %s. Access from thread %s not allowed.",
				thread.getName(), currentThread.getName()
			));
		}
	}
}
