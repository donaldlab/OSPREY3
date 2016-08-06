package edu.duke.cs.osprey.gpu;

import java.util.Set;

import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLDevice;

public class Gpu {
	
	private CLDevice device;
	
	public Gpu(CLDevice device) {
		this.device = device;
	}
	
	public CLDevice getDevice() {
		return device;
	}
	
	public CLCommandQueue makeQueue() {
		return makeQueue(false);
	}
	
	public CLCommandQueue makeQueue(boolean useProfiling) {
		if (useProfiling) {
			return device.createCommandQueue(CLCommandQueue.Mode.PROFILING_MODE);
		} else {
			return device.createCommandQueue();
		}
	}
	
	public boolean supportsDoubles() {
		Set<String> extensions = device.getExtensions();
		return extensions.contains("cl_khr_fp64") || extensions.contains("cl_amd_fp64");
	}
	
	@Override
	public String toString() {
		return device.getName();
	}
}
