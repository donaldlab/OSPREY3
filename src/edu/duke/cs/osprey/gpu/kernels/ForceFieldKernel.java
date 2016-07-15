package edu.duke.cs.osprey.gpu.kernels;

import java.io.IOException;
import java.nio.DoubleBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLMemory;

import edu.duke.cs.osprey.gpu.BoundKernel;
import edu.duke.cs.osprey.gpu.Gpu;
import edu.duke.cs.osprey.gpu.Kernel;

public class ForceFieldKernel extends Kernel<ForceFieldKernel.Bound> {
	
	public ForceFieldKernel()
	throws IOException {
		super("forcefield.cl", "calc");
	}
	
	@Override
	public Bound bind(Gpu gpu) {
		return new Bound(this, gpu);
	}
	
	public static class Bound extends BoundKernel<Bound> {
		
		private CLBuffer<DoubleBuffer> bufIn;
		private CLBuffer<DoubleBuffer> bufOut;
		
		public Bound(Kernel<ForceFieldKernel.Bound> kernel, Gpu gpu) {
			super(kernel, gpu);
		}
		
		@Override
		protected void initBuffers(int workSize) {
			bufIn = makeOrIncreaseBuffer(bufIn, workSize, CLMemory.Mem.READ_ONLY);
			bufOut = makeOrIncreaseBuffer(bufOut, workSize, CLMemory.Mem.WRITE_ONLY);
			getKernel().getCLKernel().putArgs(bufIn, bufOut);
		}
		
		public DoubleBuffer getIn() {
			return bufIn.getBuffer();
		}
		
		public DoubleBuffer getOut() {
			return bufOut.getBuffer();
		}

		@Override
		public void uploadAsync() {
			uploadBufferAsync(bufIn);
		}

		@Override
		public void downloadSync() {
			downloadBufferSync(bufOut);
		}
	}
}
