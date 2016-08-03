package edu.duke.cs.osprey.gpu;

import java.util.Set;

import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLPlatform;

public class Diagnostics {
	
	public static void main(String[] args) {
		
		// TODO: sometimes you have to `sudo clinfo` to "wake up" the OpenCL implementation
		// before it will work for non-root users. =(
		// need to find a workaround for that...
		
		// show platforms
		System.out.println("All OpenCL platforms:");
		for (int i=0; i<CLPlatform.listCLPlatforms().length; i++) {
			CLPlatform platform = CLPlatform.listCLPlatforms()[i];
			System.out.println("\n" + platform.getName());
			System.out.println("\tVendor:       " + platform.getVendor());
			System.out.println("\tVersion:      " + platform.getVersion());
			
			// show gpus
			System.out.println("\tDevices:");
			for (CLDevice device : platform.listCLDevices()) {
				System.out.println("\n\t" + device.getName());
				System.out.println("\t\tprocessors (not cores):  " + device.getMaxComputeUnits());
				System.out.println("\t\tclock speed:             " + device.getMaxClockFrequency() + " Mhz");
				System.out.println("\t\tglobal memory:           " + device.getGlobalMemSize()/1024/1024 + " MiB");
				System.out.println("\t\tlocal memory:            " + device.getLocalMemSize()/1024 + " KiB");
				System.out.println("\t\tsamplers:                " + device.getMaxSamplers());
				System.out.println("\t\t64-bit double support:   " + supportsDoubles(device));
				System.out.println("\t\tMax workgroup size:      " + device.getMaxWorkGroupSize());
			}
		}
		
		// show the usable GPUs
		System.out.println("\nGPUs usable by Osprey:");
		Gpus gpus = Gpus.get();
		for (Gpu gpu : gpus.getGpus()) {
			System.out.print("\t" + gpu);
			if (gpus.getGpus().size() > 1 && gpu == gpus.getBestGpu()) {
				System.out.print("   (best GPU)");
			}
			System.out.println();
		}
	}
	
	public static boolean supportsDoubles(CLDevice device) {
		Set<String> extensions = device.getExtensions();
		return extensions.contains("cl_khr_fp64") || extensions.contains("cl_amd_fp64");
	}
}
