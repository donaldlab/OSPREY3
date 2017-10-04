package edu.duke.cs.osprey.gpu.cuda;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jcuda.driver.CUdevice;
import jcuda.driver.JCudaDriver;

public class Gpus {
	
	public static boolean printGpuInfo = true;
	
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
		
		print("Discovering CUDA GPUs...");
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
		} catch (UnsatisfiedLinkError ex) {
			StringWriter buf = new StringWriter();
			ex.printStackTrace(new PrintWriter(buf));
			print(buf.toString());
		} finally {
			if (gpus.isEmpty()) {
				print(" none found\n");
			} else {
				print(" found " + gpus.size() + "\n");
			}
		}
	}
	
	public List<Gpu> getGpus() {
		return Collections.unmodifiableList(gpus);
	}
	
	private void print(String msg) {
		if (printGpuInfo) {
			System.out.print(msg);
		}
	}
}
