package edu.duke.cs.osprey.gpu.cuda;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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
		
		// according to docs, init flags must always be zero
		JCudaDriver.setExceptionsEnabled(true);
		JCudaDriver.cuInit(0);
		
		// how many gpus are there?
		int[] ints = new int[1];
		JCudaDriver.cuDeviceGetCount(ints);
		int count = ints[0];
		
		// get the ones that have double support
		gpus = new ArrayList<>();
		for (int i=0; i<count; i++) {
			
			CUdevice device = new CUdevice();
			JCudaDriver.cuDeviceGet(device, i);
			Gpu gpu = new Gpu(device);
			
			if (gpu.supportsDoubles()) {
				gpus.add(gpu);
			}
		}
		
		if (gpus.isEmpty()) {
			System.out.println(" none found");
		} else {
			System.out.println(" found " + gpus.size());
		}
	}
	
	public List<Gpu> getGpus() {
		
		// AAO 2017: first sort gpus by decreasing order of free memory as a
		// crude load-balancing strategy. otherwise, gpus are always hit in
		// numerical device order for multiple osprey instances, leading to 
		// unnecessary resource contention.
		Collections.sort(gpus, new Comparator<Gpu>() {
			@Override
			public int compare(Gpu o1, Gpu o2) {
				return o1.getFreeMemory() > o2.getFreeMemory() ? -1 : 1;
			}
		});
		
		return Collections.unmodifiableList(gpus);
	}
}
