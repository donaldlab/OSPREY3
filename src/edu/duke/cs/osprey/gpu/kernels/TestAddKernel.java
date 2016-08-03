package edu.duke.cs.osprey.gpu.kernels;

import java.io.IOException;
import java.nio.DoubleBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLMemory;

import edu.duke.cs.osprey.gpu.BoundKernel;
import edu.duke.cs.osprey.gpu.Gpu;
import edu.duke.cs.osprey.gpu.Kernel;

public class TestAddKernel extends Kernel<TestAddKernel.Bound> {
	
	public TestAddKernel()
	throws IOException {
		super("test.cl", "add");
	}
	
	@Override
	public Bound bind(Gpu gpu, boolean useProfiling) {
		return new Bound(this, gpu, useProfiling);
	}
	
	public static class Bound extends BoundKernel<Bound> {
		
		private CLBuffer<DoubleBuffer> bufA;
		private CLBuffer<DoubleBuffer> bufB;
		private CLBuffer<DoubleBuffer> bufOut;
		
		private int workSize;
		private int groupSize;
		
		public Bound(Kernel<TestAddKernel.Bound> kernel, Gpu gpu, boolean useProfiling) {
			super(kernel, gpu, useProfiling);
		}
		
		public DoubleBuffer getA() {
			return bufA.getBuffer();
		}
		
		public DoubleBuffer getB() {
			return bufB.getBuffer();
		}
		
		public DoubleBuffer getOut() {
			return bufOut.getBuffer();
		}
		
		public void setArgs(int workSize) {
			groupSize = getMaxGroupSize();
			this.workSize = roundUpWorkSize(workSize, groupSize);
			bufA = makeOrIncreaseBuffer(bufA, workSize, CLMemory.Mem.READ_ONLY);
			bufB = makeOrIncreaseBuffer(bufB, workSize, CLMemory.Mem.READ_ONLY);
			bufOut = makeOrIncreaseBuffer(bufOut, workSize, CLMemory.Mem.WRITE_ONLY);
		}
		
		public void runAsync() {
			getKernel().getCLKernel()
				.putArg(bufA)
				.putArg(bufB)
				.putArg(bufOut)
				.rewind();
			runAsync(workSize, groupSize);
		}
		
		public void uploadAsync() {
			uploadBufferAsync(bufA);
			uploadBufferAsync(bufB);
		}

		public void downloadSync() {
			downloadBufferSync(bufOut);
		}
	}
}
