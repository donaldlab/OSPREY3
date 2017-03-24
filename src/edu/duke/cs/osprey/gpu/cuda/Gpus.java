package edu.duke.cs.osprey.gpu.cuda;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jcuda.driver.CUdevice;
import jcuda.driver.JCudaDriver;

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
		
		System.out.print("Discovering CUDA GPUs...");
		gpus = new ArrayList<>();

		try {
			// according to docs, init flags must always be zero
			JCudaDriver.setExceptionsEnabled(true);
			JCudaDriver.cuInit(0);
			
			// how many gpus are there?
			int[] ints = new int[1];
			JCudaDriver.cuDeviceGetCount(ints);
			int count = ints[0];
			
			// get the ones that have double support
			for (int i=0; i<count; i++) {
				
				CUdevice device = new CUdevice();
				JCudaDriver.cuDeviceGet(device, i);
				Gpu gpu = new Gpu(device);
				
				if (gpu.supportsDoubles()) {
					gpus.add(gpu);
				}
			}
		} finally {
			if (gpus.isEmpty()) {
				System.out.println(" none found");
			} else {
				System.out.println(" found " + gpus.size());
			}
		}
	}
	
	public List<Gpu> getGpus() {
		return Collections.unmodifiableList(gpus);
	}
}
