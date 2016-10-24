package edu.duke.cs.osprey.gpu.cuda;

import java.util.Arrays;
import java.util.List;

public class ContextPool {
	
	private List<Gpu> gpus;
	private boolean[] checkedOut;
	
	public ContextPool() {
		this(Gpus.get().getGpus().size());
	}
	
	public ContextPool(int numGpus) {
		
		// make sure we don't try to use too many gpus
		List<Gpu> gpus = Gpus.get().getGpus();
		numGpus = Math.min(numGpus, gpus.size());
		this.gpus = gpus.subList(0, numGpus);
		
		// initially, all gpus are available
		checkedOut = new boolean[numGpus];
		Arrays.fill(checkedOut, false);
		
		System.out.println(String.format("ContextPool: using %d gpu(s)", numGpus));
	}
	
	public synchronized Context checkout() {
		
		// NOTE: can't cache contexts even if the same thread is asking over and over
		// there's no way to clean them up from another thread
		// we *need* to cleanup the contexts at release() time
		// unless we add an explicit per-thread cleanup method
		
		// find an available gpu
		for (int i=0; i<gpus.size(); i++) {
			if (!checkedOut[i]) {
				checkedOut[i] = true;
				return new Context(gpus.get(i));
			}
		}
		
		throw new IllegalStateException(String.format("no more queues to checkout, %d already used", gpus.size()));
	}
	
	public synchronized void release(Context context) {
		
		context.cleanup();
		
		for (int i=0; i<gpus.size(); i++) {
			if (context.getGpu() == gpus.get(i)) {
				checkedOut[i] = false;
			}
		}
	}
}
