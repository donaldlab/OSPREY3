package edu.duke.cs.osprey.gpu.opencl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLException;
import com.jogamp.opencl.CLPlatform;

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
		
		System.out.print("Discovering OpenCL GPUs...");
		
		// get the gpus that support doubles
		gpus = new ArrayList<>();
		for (int i=0; i<CLPlatform.listCLPlatforms().length; i++) {
			CLPlatform platform = CLPlatform.listCLPlatforms()[i];
			for (CLDevice device : platform.listCLDevices()) {
				
				// CPUs are slow, don't even bother
				if (device.getType() == CLDevice.Type.CPU) {
					continue;
				}

				try {
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
				} catch (CLException.CLInvalidDeviceException ex) {
					// OpenCL doesn't work on this device, skip it
				}
			}
		}
		
		if (gpus.isEmpty()) {
			System.out.println(" none found");
		} else {
			System.out.println(" found " + gpus.size());
		}
	}
	
	public List<Gpu> getGpus() {
		return Collections.unmodifiableList(gpus);
	}
}
