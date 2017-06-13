package edu.duke.cs.osprey.gpu.cuda;

import java.nio.Buffer;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.tpie.Cleaner;
import edu.duke.cs.tpie.Cleaner.GarbageDetectable;

public class BufferPool<T extends Buffer> implements GarbageDetectable {
	
	public static interface BufferFactory<T extends Buffer> extends Factory<CUBuffer<T>,Integer> {
		// nothing else to do
	}
	
	public static interface BufferExpander<T extends Buffer> {
		CUBuffer<T> expand(CUBuffer<T> buf, int size);
	}

	public final BufferFactory<T> factory;
	public final BufferExpander<T> expand;
	
	private TreeMap<Integer,List<CUBuffer<T>>> buffers;
	private Cleaner.Cleanable cleaner;
	
	public BufferPool(BufferFactory<T> factory, BufferExpander<T> expand) {
		
		this.factory = factory;
		this.expand = expand;
		this.buffers = new TreeMap<>();
		
		// setup the cleaner
		// NOTE: copy the buffers ref, so the cleaner doesn't hold a strong reference to this
		final TreeMap<Integer,List<CUBuffer<T>>> buffers = this.buffers;
		cleaner = () -> {
			for (List<CUBuffer<T>> bufs : buffers.values()) {
				for (CUBuffer<T> buf : bufs) {
					buf.cleanup();
				}
			}
			buffers.clear();
		};
		Cleaner.addCleaner(this, cleaner);
	}
	
	public void cleanup() {
		cleaner.clean();
	}
	
	private List<CUBuffer<T>> getBuffers(int size) {
		List<CUBuffer<T>> bufs = buffers.get(size);
		if (bufs == null) {
			bufs = new ArrayList<>();
			buffers.put(size, bufs);
		}
		return bufs;
	}
	
	private void add(CUBuffer<T> buf) {
		getBuffers(buf.size()).add(buf);
	}
	
	private CUBuffer<T> remove(int size) {
		
		if (!buffers.containsKey(size)) {
			return null;
		}
		
		List<CUBuffer<T>> bufs = getBuffers(size);
		if (bufs.isEmpty()) {
			
			// this shouldn't happen, but just in case...
			buffers.remove(size);
			return null;
		}
			
		CUBuffer<T> buf = bufs.remove(bufs.size() - 1);
		
		if (bufs.isEmpty()) {
			buffers.remove(size);
		}
		
		return buf;
	}
	
	public CUBuffer<T> checkout(int size) {
		
		// look for a buffer that's big enough
		Integer sizeOrBigger = buffers.ceilingKey(size);
		if (sizeOrBigger != null) {
			return prep(remove(sizeOrBigger), size);
		}
		
		// look for a buffer that's too small
		Integer sizeOrSmaller = buffers.floorKey(size);
		if (sizeOrSmaller != null) {
			
			// this buffer is too small, so expand it
			// that way, we don't have a bunch of tiny buffers cluttering up memory
			CUBuffer<T> buf = remove(sizeOrSmaller);
			buf = expand.expand(buf, size);
			return prep(buf, size);
		}
		
		// make a new one
		return prep(factory.make(size), size);
	}
	
	public void release(CUBuffer<T> buf) {
		add(buf);
	}
	
	private CUBuffer<T> prep(CUBuffer<T> buf, int size) {
		
		// make the buffer ready to use
		buf.getHostBuffer().clear();
		
		// just in case
		assert (buf.getHostBuffer().capacity() >= size);
		
		return buf;
	}
}
