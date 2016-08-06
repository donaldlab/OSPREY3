package edu.duke.cs.osprey.gpu;

import java.util.ArrayList;
import java.util.List;

import com.jogamp.opencl.CLCommandQueue;

public class GpuQueuePool {
	
	private int queuesPerGpu;
	private List<List<CLCommandQueue>> queues;
	private int nextIndex;
	
	public GpuQueuePool(int numGpus, int queuesPerGpu) {
		this(numGpus, queuesPerGpu, false);
	}
	
	public GpuQueuePool(int numGpus, int queuesPerGpu, boolean useProfiling) {
		
		this.queuesPerGpu = queuesPerGpu;
		
		// make sure we don't try to use too many gpus
		List<Gpu> gpus = Gpus.get().getGpus();
		numGpus = Math.min(numGpus, gpus.size());
		
		// make n queues for each gpu
		queues = new ArrayList<>(numGpus);
		for (int i=0; i<numGpus; i++) {
			Gpu gpu = gpus.get(i);
			List<CLCommandQueue> queuesAtGpu = new ArrayList<>();
			for (int j=0; j<queuesPerGpu; j++) {
				queuesAtGpu.add(gpu.makeQueue(useProfiling));
			}
			queues.add(queuesAtGpu);
		}
		
		nextIndex = 0;
	}
	
	public int getNumGpus() {
		return queues.size();
	}
	
	public int getQueuesPerGpu() {
		return queuesPerGpu;
	}
	
	public int getNumQueues() {
		return queues.size()*queuesPerGpu;
	}
	
	public CLCommandQueue getQueue(int index) {
		
		// visit the first queue on each gpu first, then the second, etc...
		int gpuIndex = index%queues.size();
		int queueIndex = index/queues.size();
		
		return getQueue(gpuIndex, queueIndex);
	}
	
	public CLCommandQueue getQueue(int gpuIndex, int queueIndex) {
		return queues.get(gpuIndex).get(queueIndex);
	}
	
	public CLCommandQueue getRoundRobinQueue() {
		return getQueue(nextIndex++ % getNumQueues());
	}
}
