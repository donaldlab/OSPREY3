/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.gpu.cuda;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import edu.duke.cs.osprey.tools.FileTools.ResourcePathRoot;
import jcuda.Pointer;
import jcuda.driver.CUcontext;
import jcuda.driver.CUctx_flags;
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
		//int flags = CUctx_flags.CU_CTX_SCHED_YIELD;
		//int flags = CUctx_flags.CU_CTX_SCHED_SPIN;
		int flags = CUctx_flags.CU_CTX_SCHED_BLOCKING_SYNC;
		JCudaDriver.cuCtxCreate(context, flags, gpu.getDevice());
		
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
			File kernelFile = new ResourcePathRoot("/gpuKernels/cuda").extractToTempFile(String.format("%s.bin", name));
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
	
	public void uploadAsync(CUdeviceptr pdBuf, Pointer phBuf, long numBytes, GpuStream stream) {
		JCudaDriver.cuMemcpyHtoDAsync(pdBuf, phBuf, numBytes, stream.getStream());
	}
	
	public void downloadAsync(Pointer phBuf, CUdeviceptr pdBuf, long numBytes, GpuStream stream) {
		JCudaDriver.cuMemcpyDtoHAsync(phBuf, pdBuf, numBytes, stream.getStream());
	}
	
	public void pinBuffer(Pointer phBuf, long numBytes) {
		if (numBytes <= 0) {
			throw new IllegalArgumentException("bad buffer size: " + numBytes + " bytes");
		}
		JCudaDriver.cuMemHostRegister(phBuf, numBytes, 0);
	}
	
	public void unpinBuffer(Pointer phBuf) {
		JCudaDriver.cuMemHostUnregister(phBuf);
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
		JCudaDriver.cuCtxSetCurrent(context);
	}
	
	public synchronized void cleanup() {
		try {
			for (CUmodule kernel : kernels.values()) {
				JCudaDriver.cuModuleUnload(kernel);
			}
			kernels.clear();
			
			JCudaDriver.cuCtxDestroy(context);
		} catch (Throwable t) {
			t.printStackTrace(System.err);
		}
	}
}
