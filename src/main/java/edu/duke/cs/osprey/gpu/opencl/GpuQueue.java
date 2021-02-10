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

import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;

public class GpuQueue {
	
	private Gpu gpu;
	private CLCommandQueue queue;
	private CLContext separateContext;
	
	public GpuQueue(Gpu gpu) {
		this(gpu, false, false);
	}
	
	public GpuQueue(Gpu gpu, boolean useProfiling, boolean makeSeparateContext) {
		
		if (makeSeparateContext) {
			separateContext = CLContext.create(gpu.getDevice());
			this.gpu = new Gpu(separateContext.getDevices()[0]);
		} else {
			separateContext = null;
			this.gpu = gpu;
		}
		
		CLDevice device = this.gpu.getDevice();
		if (useProfiling) {
			queue = device.createCommandQueue(CLCommandQueue.Mode.PROFILING_MODE);
		} else {
			queue = device.createCommandQueue();
		}
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	public CLCommandQueue getCLQueue() {
		return queue;
	}
	
	public void cleanup() {
		
		if (queue != null) {
			queue.release();
			queue = null;
		}
		
		if (separateContext != null) {
			separateContext.release();
			separateContext = null;
		}
	}
	
	@Override
	protected void finalize()
	throws Throwable {
		try {
			if (queue != null || separateContext != null) {
				System.err.println("WARNING: " + getClass().getName() + " was garbage collected, but not cleaned up. Attempting cleanup now");
				cleanup();
			}
		} finally {
			super.finalize();
		}
	}

	public boolean isProfilingEnabled() {
		return queue.isProfilingEnabled();
	}
	
	public void flush() {
		queue.flush();
	}

	public void waitForGpu() {
		queue.finish();
	}
}
