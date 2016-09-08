package edu.duke.cs.osprey.gpu;

import java.util.Set;

import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLPlatform;

public class Diagnostics {
	
	public static void main(String[] args) {
		
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
				System.out.println(String.format("\t\t%-30s %b", "64-bit double support:", supportsDoubles(device)));
				System.out.println(String.format("\t\t%-30s %d", "compute units:", device.getMaxComputeUnits()));
				System.out.println(String.format("\t\t%-30s %d", "work items per compute unit:", device.getMaxWorkGroupSize()/device.getMaxComputeUnits()));
				System.out.println(String.format("\t\t%-30s %d", "max workgroup size:", device.getMaxWorkGroupSize()));
				System.out.println(String.format("\t\t%-30s %d Mhz", "clock speed:", device.getMaxClockFrequency()));
				System.out.println(String.format("\t\t%-30s %d MiB", "global memory:", device.getGlobalMemSize()/1024/1024));
				System.out.println(String.format("\t\t%-30s %d KiB", "local memory:", device.getLocalMemSize()/1024));
				System.out.println(String.format("\t\t%-30s %d MiB", "max allocation:", device.getMaxMemAllocSize()/1024/1024));
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
