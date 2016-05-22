package edu.duke.cs.osprey.gpu;

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
}
