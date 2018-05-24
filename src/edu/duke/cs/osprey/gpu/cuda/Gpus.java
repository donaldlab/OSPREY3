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
