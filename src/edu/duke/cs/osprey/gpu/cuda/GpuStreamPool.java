package edu.duke.cs.osprey.gpu.cuda;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class GpuStreamPool {
	
	private int numStreamsPerGpu;
	private List<Context> contexts;
	private List<List<GpuStream>> streamsByGpu;
	private List<GpuStream> streams;
	private boolean[] checkedOut;
		
	public GpuStreamPool() {
		this(1);
	}
	
	public GpuStreamPool(int queuesPerGpu) {
		this(Gpus.get().getGpus().size(), queuesPerGpu);
	}
	
	public GpuStreamPool(int numGpus, int queuesPerGpu) {
		this(numGpus, queuesPerGpu, false);
	}
	
	public GpuStreamPool(int numGpus, int numStreamsPerGpu, boolean useProfiling) {
		
		this.numStreamsPerGpu = numStreamsPerGpu;
		
		// make sure we don't try to use too many gpus
		List<Gpu> gpus = new ArrayList<Gpu>(Gpus.get().getGpus());

		// AAO 2017: first sort gpus by decreasing order of free memory as a
		// crude load-balancing strategy. otherwise, gpus are always hit in
		// numerical device order for multiple osprey instances, leading to 
		// unnecessary resource contention.
		Collections.sort(gpus, new Comparator<Gpu>() {
			@Override
			public int compare(Gpu o1, Gpu o2) {
				return o1.getFreeMemory() > o2.getFreeMemory() ? -1 : 1;
			}
		});

		numGpus = Math.min(numGpus, gpus.size());
		
		// make contexts for all the gpus
		contexts = new ArrayList<>();
		for (Gpu gpu : gpus) {
			Context context = new Context(gpu);
			contexts.add(context);
		}
		
		// make n stream for each gpu
		streamsByGpu = new ArrayList<>(numGpus);
		for (int i=0; i<numGpus; i++) {
			Context context = contexts.get(i);
			context.attachCurrentThread();
			List<GpuStream> queuesAtGpu = new ArrayList<>();
			for (int j=0; j<numStreamsPerGpu; j++) {
				queuesAtGpu.add(new GpuStream(context));
			}
			streamsByGpu.add(queuesAtGpu);
		}
		
		// flatten the streams into a list
		// visit the first stream on each gpu first, then the second, etc...
		streams = new ArrayList<>(numGpus*numStreamsPerGpu);
		for (int i=0; i<numStreamsPerGpu; i++) {
			for (int j=0; j<numGpus; j++) {
				streams.add(streamsByGpu.get(j).get(i));
			}
		}
		
		// initially, all streams are available
		checkedOut = new boolean[streams.size()];
		Arrays.fill(checkedOut, false);
		
		System.out.println(String.format("GpuStreamPool: using %d stream(s) across %d gpu(s)", streams.size(), numGpus));
	}
	
	public int getNumGpus() {
		return streamsByGpu.size();
	}
	
	public int getNumStreamsPerGpu() {
		return numStreamsPerGpu;
	}
	
	public int getNumStreams() {
		return streams.size();
	}
	
	public synchronized GpuStream checkout() {
		
		// find an available queue
		for (int i=0; i<streams.size(); i++) {
			if (!checkedOut[i]) {
				checkedOut[i] = true;
				GpuStream stream = streams.get(i);
				stream.getContext().attachCurrentThread();
				return stream;
			}
		}
		
		throw new IllegalStateException(String.format("no more streams to checkout, %d already used", streams.size()));
	}
	
	public synchronized void release(GpuStream stream) {
		
		for (int i=0; i<streams.size(); i++) {
			if (streams.get(i) == stream) {
				checkedOut[i] = false;
			}
		}
	}

	public void cleanup() {
		for (GpuStream stream : streams) {
			stream.cleanup();
		}
		streams.clear();
		for (Context context : contexts) {
			context.cleanup();
		}
		contexts.clear();
	}
}
