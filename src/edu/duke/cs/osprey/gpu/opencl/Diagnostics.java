/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.gpu.opencl;

import java.util.Set;

import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLException.CLPlatformNotFoundKhrException;
import com.jogamp.opencl.CLPlatform;

public class Diagnostics {
	
	public static void main(String[] args) {
		
		try {
		
			// show platforms
			System.out.println("Looking for OpenCL GPUs...");
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
			System.out.println("\nOpenCL GPUs usable by Osprey:");
			Gpus gpus = Gpus.get();
			for (Gpu gpu : gpus.getGpus()) {
				System.out.println("\t" + gpu);
			}
			
		} catch (UnsatisfiedLinkError ex) {
			System.out.println("OpenCL does not appear to be installed. No GPUs found.");
		} catch (CLPlatformNotFoundKhrException ex) {
			System.out.println("No OpenCL GPUs found.");
		}
	}
	
	public static boolean supportsDoubles(CLDevice device) {
		Set<String> extensions = device.getExtensions();
		return extensions.contains("cl_khr_fp64") || extensions.contains("cl_amd_fp64");
	}
}
