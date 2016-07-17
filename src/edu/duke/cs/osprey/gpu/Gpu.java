package edu.duke.cs.osprey.gpu;

import java.util.Set;

import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLDevice;

public class Gpu {
	
	private CLDevice device;
	private CLCommandQueue queue;
	
	public Gpu(CLDevice device) {
		this.device = device;
		this.queue = device.createCommandQueue();
	}
	
	public CLDevice getDevice() {
		return device;
	}
	
	public CLCommandQueue getQueue() {
		return queue;
	}
	
	public boolean supportsDoubles() {
		Set<String> extensions = device.getExtensions();
		return extensions.contains("cl_khr_fp64") || extensions.contains("cl_amd_fp64");
	}
	
	@Override
	public String toString() {
		return device.getName();
	}
	
	public void cleanup() {
		queue.release();
	}
}
