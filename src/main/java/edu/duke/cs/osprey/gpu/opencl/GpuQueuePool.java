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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GpuQueuePool {
	
	private int numQueuesPerGpu;
	private List<List<GpuQueue>> queuesByGpu;
	private List<GpuQueue> queues;
	private boolean[] checkedOut;
	
	public GpuQueuePool() {
		this(1);
	}
	
	public GpuQueuePool(int queuesPerGpu) {
		this(Gpus.get().getGpus().size(), queuesPerGpu);
	}
	
	public GpuQueuePool(int numGpus, int queuesPerGpu) {
		this(numGpus, queuesPerGpu, false);
	}
	
	public GpuQueuePool(int numGpus, int numQueuesPerGpu, boolean useProfiling) {
		
		this.numQueuesPerGpu = numQueuesPerGpu;
		
		// make sure we don't try to use too many gpus
		List<Gpu> gpus = Gpus.get().getGpus();
		numGpus = Math.min(numGpus, gpus.size());
		
		// make n queues for each gpu
		queuesByGpu = new ArrayList<>(numGpus);
		for (int i=0; i<numGpus; i++) {
			Gpu gpu = gpus.get(i);
			List<GpuQueue> queuesAtGpu = new ArrayList<>();
			for (int j=0; j<numQueuesPerGpu; j++) {
				queuesAtGpu.add(new GpuQueue(gpu, useProfiling, false));
			}
			queuesByGpu.add(queuesAtGpu);
		}
		
		// flatten the queues into a list
		// visit the first queue on each gpu first, then the second, etc...
		queues = new ArrayList<>(numGpus*numQueuesPerGpu);
		for (int i=0; i<numQueuesPerGpu; i++) {
			for (int j=0; j<numGpus; j++) {
				queues.add(queuesByGpu.get(j).get(i));
			}
		}
		
		// initially, all queues are available
		checkedOut = new boolean[queues.size()];
		Arrays.fill(checkedOut, false);
		
		System.out.println(String.format("GpuQueuePool: using %d command queue(s) across %d gpu(s)", queues.size(), numGpus));
	}
	
	public int getNumGpus() {
		return queuesByGpu.size();
	}
	
	public int getNumQueuesPerGpu() {
		return numQueuesPerGpu;
	}
	
	public int getNumQueues() {
		return queues.size();
	}
	
	public synchronized GpuQueue checkout() {
		
		// find an available queue
		for (int i=0; i<queues.size(); i++) {
			if (!checkedOut[i]) {
				checkedOut[i] = true;
				return queues.get(i);
			}
		}
		
		throw new IllegalStateException(String.format("no more queues to checkout, %d already used", queues.size()));
	}
	
	public synchronized void release(GpuQueue queue) {
		
		// if we're already cleaned up, cleanup the queue now
		if (queues == null) {
			queue.cleanup();
			return;
		}
		
		for (int i=0; i<queues.size(); i++) {
			if (queues.get(i) == queue) {
				checkedOut[i] = false;
			}
		}
	}

	public void cleanup() {
		if (queues != null) {
			for (GpuQueue queue : queues) {
				queue.cleanup();
			}
			queues = null;
			queuesByGpu = null;
		}
	}
	
	@Override
	protected void finalize()
	throws Throwable {
		try {
			if (queues != null) {
				System.err.println("WARNING: " + getClass().getName() + " was garbage collected, but not cleaned up. Attempting cleanup now");
				cleanup();
			}
		} finally {
			super.finalize();
		}
	}
}
