package edu.duke.cs.osprey.gpu;

import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLPlatform;

public class MainDiagnostics {
	
	public static void main(String[] args) {
		
		// show platforms
		System.out.println("platforms:");
		for (CLPlatform platform : CLPlatform.listCLPlatforms()) {
			System.out.println("\t" + platform);
			
			// show gpus
			for (CLDevice device : platform.listCLDevices()) {
				System.out.println("\t\t" + device);
				System.out.println("\t\t\tprocessors (not cores):  " + device.getMaxComputeUnits());
				System.out.println("\t\t\tclock speed:             " + device.getMaxClockFrequency() + " Mhz");
				System.out.println("\t\t\tmemory:                  " + device.getGlobalMemSize()/1024/1024 + " MiB");
				System.out.println("\t\t\tsamplers:                " + device.getMaxSamplers());
			}
		}
	}
}
