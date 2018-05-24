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

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.Random;

import org.junit.Test;

import edu.duke.cs.osprey.gpu.opencl.GpuQueue;
import edu.duke.cs.osprey.gpu.opencl.Gpus;
import edu.duke.cs.osprey.gpu.opencl.kernels.TestFancyKernel;

public class TestGpu {
	
	@Test
	public void testAccuracy()
	throws Exception {
		
		// init random doubles
		final int n = 1024 * 1024;
		double[] a = new double[n];
		double[] b = new double[n];
		double[] out = new double[n];
		Random rand = new Random(12345);
		for (int i=0; i<n; i++) {
			a[i] = rand.nextDouble();
			b[i] = rand.nextDouble();
			out[i] = Math.sqrt(a[i]*a[i] + b[i]*b[i]);
		}
		
		final int NumRuns = 10;
		
		GpuQueue queue = new GpuQueue(Gpus.get().getGpus().get(0));
		TestFancyKernel kernel = new TestFancyKernel(queue);
		
		// copy data to buffers
		kernel.setArgs(n);
		for (int i=0; i<n; i++) {
			kernel.getA().put(a[i]);
			kernel.getB().put(b[i]);
		}
		kernel.getA().rewind();
		kernel.getB().rewind();
		
		// upload args to gpu
		kernel.uploadSync();
		
		for (int i=0; i<NumRuns; i++) {
			
			// run the kernel
			kernel.runAsync();
			kernel.downloadSync();
			
			// check the results for accuracy
			for (int j=0; j<n; j++) {
				double gpuVal = kernel.getOut().get(j);
				double err = Math.abs(out[j] - gpuVal);
				assertThat(err, lessThanOrEqualTo(1e-15));
			}
		}
		
		kernel.cleanup();
		queue.cleanup();
	}
}
