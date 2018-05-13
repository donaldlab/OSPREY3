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
