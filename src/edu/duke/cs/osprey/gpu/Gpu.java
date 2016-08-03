package edu.duke.cs.osprey.gpu;

import java.util.Set;

import com.jogamp.opencl.CLDevice;

public class Gpu {
	
	private CLDevice device;
	
	public Gpu(CLDevice device) {
		this(device, false);
	}
	
	public Gpu(CLDevice device, boolean useProfiling) {
		this.device = device;
	}
	
	public CLDevice getDevice() {
		return device;
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
