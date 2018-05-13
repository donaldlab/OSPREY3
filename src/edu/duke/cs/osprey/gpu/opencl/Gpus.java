/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

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
