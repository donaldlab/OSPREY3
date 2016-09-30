package edu.duke.cs.osprey.gpu;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;

public class Gpus {
	
	private static Gpus instance;
	
	static {
		instance = null;
	}
	
	public static Gpus get() {
		if (instance == null) {
			instance = new Gpus();
		}
		return instance;
	}
	
	private List<Gpu> gpus;
	
	private Gpus() {
		
		System.out.println("Discovering GPUs...");
		
		// get the gpus that support doubles
		CLContext context = CLContext.create();
		gpus = new ArrayList<>();
		for (CLDevice device : context.getDevices()) {
			if (device.getType() == CLDevice.Type.GPU) {
				
				// make an independent context for each gpu
				CLContext gpuContext = CLContext.create(device);
				CLDevice isolatedDevice = gpuContext.getDevices()[0];
				
				Gpu gpu = new Gpu(isolatedDevice);
				if (gpu.supportsDoubles()) {
					
					// we can use this gpu! save it and its context
					gpus.add(gpu);
					
				} else {
					
					// can't use this gpu, clean it up
					gpuContext.release();
				}
			}
		}
	}
	
	public List<Gpu> getGpus() {
		return Collections.unmodifiableList(gpus);
	}
	
	public Gpu getBestGpu() {
		
		Gpu bestGpu = null;
		long bestScore = 0;
		
		for (Gpu gpu : gpus) {
			long score = gpu.getDevice().getMaxComputeUnits() * gpu.getDevice().getMaxClockFrequency();
			if (score > bestScore) {
				bestScore = score;
				bestGpu = gpu;
			}
		}
		
		return bestGpu;
	}
}
