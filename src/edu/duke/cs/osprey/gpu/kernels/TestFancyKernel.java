package edu.duke.cs.osprey.gpu.kernels;

import java.io.IOException;
import java.nio.DoubleBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLMemory;

import edu.duke.cs.osprey.gpu.BoundKernel;
import edu.duke.cs.osprey.gpu.Gpu;
import edu.duke.cs.osprey.gpu.Kernel;

public class TestFancyKernel extends Kernel<TestFancyKernel.Bound> {
	
	public TestFancyKernel()
	throws IOException {
		super("test.cl", "fancy");
	}
	
	@Override
	public Bound bind(Gpu gpu) {
		return new Bound(this, gpu);
	}
	
	public static class Bound extends BoundKernel<Bound> {
		
		private CLBuffer<DoubleBuffer> bufA;
		private CLBuffer<DoubleBuffer> bufB;
		private CLBuffer<DoubleBuffer> bufOut;
		
		public Bound(Kernel<TestFancyKernel.Bound> kernel, Gpu gpu) {
			super(kernel, gpu);
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
			workSize = roundUpWorkSize(workSize);
			setWorkSize(workSize);
			bufA = makeOrIncreaseBuffer(bufA, workSize, CLMemory.Mem.READ_ONLY);
			bufB = makeOrIncreaseBuffer(bufB, workSize, CLMemory.Mem.READ_ONLY);
			bufOut = makeOrIncreaseBuffer(bufOut, workSize, CLMemory.Mem.WRITE_ONLY);
			getKernel().getCLKernel()
				.putArg(bufA)
				.putArg(bufB)
				.putArg(bufOut);
		}

		public void uploadSync() {
			uploadBufferAsync(bufA);
			uploadBufferAsync(bufB);
			waitForGpu();
		}

		public void downloadSync() {
			downloadBufferSync(bufOut);
		}
	}
}
