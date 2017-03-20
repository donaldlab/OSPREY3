package edu.duke.cs.osprey.gpu.cuda;

import jcuda.driver.CUcontext;
import jcuda.driver.CUdevice;
import jcuda.driver.CUdevice_attribute;
import jcuda.driver.JCudaDriver;

public class Gpu {
	
	private CUdevice device;
	private String name;
	private int[] computeVersion;
	private int warpThreads;
	private int maxBlockThreads;
	private long totalMemory;
	private long freeMemory;

	public Gpu(CUdevice device) {
		
		this.device = device;
		
		// get name
		byte[] bytes = new byte[1024];
		JCudaDriver.cuDeviceGetName(bytes, bytes.length, device);
		name = new String(bytes);
		
		// get total and free memory
		CUcontext cuCtx = new CUcontext();
		JCudaDriver.cuCtxCreate(cuCtx, 0, device);
		long[][] longs = new long[2][1];
		JCudaDriver.cuMemGetInfo(longs[0], longs[1]);
		freeMemory = longs[0][0];
		totalMemory = longs[1][0];
		JCudaDriver.cuCtxDestroy(cuCtx);
		
		// get attributes
		computeVersion = new int[] {
			getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR),
			getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR)
		};
		warpThreads = getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_WARP_SIZE);
		maxBlockThreads = getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X);
	}
	
	public CUdevice getDevice() {
		return device;
	}
	
	public int getAttribute(int attr) {
		int[] ints = new int[1];
		JCudaDriver.cuDeviceGetAttribute(ints, attr, device);
		return ints[0];
	}
	
	public String getName() {
		return name;
	}
	
	public long getFreeMemory() {
		return freeMemory;
	}
	
	public long getTotalMemory() {
		return totalMemory;
	}
	
	public int[] getComputeVersion() {
		return computeVersion;
	}
	
	public boolean isComputeVersionAtLeast(int major, int minor) {
		
		if (computeVersion[0] < major) {
			return false;
		} else if (computeVersion[0] > major) {
			return true;
		}
		
		return computeVersion[1] >= minor; 
	}
	
	public boolean supportsDoubles() {
		return isComputeVersionAtLeast(1, 3);
	}
	
	public boolean supportsDynamicParallelism() {
		return isComputeVersionAtLeast(3, 5);
	}
	
	public int getWarpThreads() {
		return warpThreads;
	}
	
	public int getMaxBlockThreads() {
		return maxBlockThreads;
	}
	
	@Override
	public String toString() {
		return getName();
	}
}
