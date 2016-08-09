package edu.duke.cs.osprey.gpu;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

import com.jogamp.opencl.CLCommandQueue;

public class GpuQueuePool {
	
	private int numQueuesPerGpu;
	private List<List<CLCommandQueue>> queuesByGpu;
	private List<CLCommandQueue> queues;
	private Deque<CLCommandQueue> availableQueues;
	
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
			List<CLCommandQueue> queuesAtGpu = new ArrayList<>();
			for (int j=0; j<numQueuesPerGpu; j++) {
				queuesAtGpu.add(gpu.makeQueue(useProfiling));
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
		availableQueues = new ArrayDeque<>(queues);
		
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
	
	public synchronized CLCommandQueue checkout() {
		
		if (availableQueues.isEmpty()) {
			throw new IllegalStateException("no more queues to checkout");
		}
		
		return availableQueues.removeFirst();
	}
	
	public synchronized void release(CLCommandQueue queue) {
		availableQueues.addLast(queue);
	}

	public void cleanup() {
		for (CLCommandQueue queue : queues) {
			queue.release();
		}
		queuesByGpu.clear();
	}
}
