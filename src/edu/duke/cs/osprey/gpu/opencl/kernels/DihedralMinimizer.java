package edu.duke.cs.osprey.gpu.opencl.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLMemory;

import edu.duke.cs.osprey.gpu.opencl.GpuQueue;
import edu.duke.cs.osprey.gpu.opencl.Kernel;

public class DihedralMinimizer {
	
	private static String SourcePath = "dihedralMinimizer.cl";
	
	public static class PoseKernel extends Kernel {
	
		private int numRotatedAtoms;
		private CLBuffer<IntBuffer> dihedralIndices;
		private CLBuffer<IntBuffer> rotatedIndices;
		private CLBuffer<ByteBuffer> args;
		
		public PoseKernel(GpuQueue queue)
		throws IOException {
			super(queue, SourcePath, "pose");
		}
		
		public void initForcefield(CLBuffer<DoubleBuffer> coords, int[] dihedralIndices, int[] rotatedIndices, CLBuffer<ByteBuffer> efuncArgs) {
			
			// allocate buffers
			CLContext context = getQueue().getCLQueue().getContext();
			this.dihedralIndices = context.createIntBuffer(dihedralIndices.length, CLMemory.Mem.READ_ONLY);
			this.rotatedIndices = context.createIntBuffer(rotatedIndices.length, CLMemory.Mem.READ_ONLY);
			
			// upload the indices
			
			this.dihedralIndices.getBuffer().rewind();
			this.dihedralIndices.getBuffer().put(dihedralIndices);
			this.dihedralIndices.getBuffer().rewind();
			uploadBufferAsync(this.dihedralIndices);
			
			this.rotatedIndices.getBuffer().rewind();
			this.rotatedIndices.getBuffer().put(rotatedIndices);
			this.rotatedIndices.getBuffer().rewind();
			uploadBufferAsync(this.rotatedIndices);
			
			// make the args buffer
			int argSize = 144;
			args = context.createByteBuffer(argSize, CLMemory.Mem.READ_WRITE);
			ByteBuffer buf = args.getBuffer();
			buf.rewind();
			for (int i=0; i<2; i++) {
				buf.putInt(0);
			}
			for (int i=0; i<16; i++) {
				buf.putDouble(0);
			}
			for (int i=0; i<8; i++) {
				buf.put((byte)0);
			}
			buf.clear(); // doesn't actually empty the buffer, just resets pos and limit
			
			// set kernel args
			getCLKernel()
				.setArg(0, coords)
				.setArg(1, this.dihedralIndices)
				.setArg(2, this.rotatedIndices)
				.setArg(3, efuncArgs)
				.setArg(4, args);
			
			numRotatedAtoms = rotatedIndices.length;
		}
		
		public void initArgs(double step, double xdmin, double xdmax, double xd, double internalSolvEnergy) {
			
			// zero out the buffer
			ByteBuffer buf = args.getBuffer();
			buf.clear();
			for (int i=0; i<buf.capacity(); i++) {
				buf.put((byte)0);
			}
			
			// set initial values
			// (see Args struct in kernel source for addresses)
			buf.putInt(0, numRotatedAtoms);
			buf.putDouble(8, xdmin);
			buf.putDouble(16, xdmax);
			buf.putDouble(24, xd);
			buf.putDouble(40, step);
			buf.putDouble(120, xd); // dihedralDegrees = xd
			buf.putDouble(128, internalSolvEnergy);
			
			// upload the args
			buf.rewind();
			uploadBufferAsync(args);
		}
		
		public CLBuffer<ByteBuffer> getArgs() {
			return args;
		}
		
		public ByteBuffer downloadArgsSync() {
			args.getBuffer().rewind();
			downloadBufferSync(args);
			return args.getBuffer();
		}
		
		public void runAsync() {
			super.runAsync(numRotatedAtoms, numRotatedAtoms);
		}
		
		@Override
		public void cleanup() {
			dihedralIndices.release();
			rotatedIndices.release();
			args.release();
		}
	}
	
	private static final int EnergyReadThreads = 32; // about a hardware workgroup/warp size seems good here
	
	public static class SearchKernel extends Kernel {
		
		public SearchKernel(GpuQueue queue, int number)
		throws IOException {
			super(queue, SourcePath, "search" + number);
		}
		
		public void init(CLBuffer<DoubleBuffer> energies, int numEnergies, CLBuffer<ByteBuffer> efuncArgs, CLBuffer<ByteBuffer> args) {
			getCLKernel()
				.setArg(0, energies)
				.setArg(1, numEnergies)
				.setArg(2, efuncArgs)
				.setArg(3, args)
				.setNullArg(4, (EnergyReadThreads + 1)*Double.BYTES);
		}
		
		public void runAsync() {
			super.runAsync(EnergyReadThreads, EnergyReadThreads);
		}
	}
	
	public static class SurfKernel extends Kernel {
		
		public SurfKernel(GpuQueue queue)
		throws IOException {
			super(queue, SourcePath, "surf");
		}
		
		public void init(CLBuffer<DoubleBuffer> energies, int numEnergies, CLBuffer<ByteBuffer> efuncArgs, CLBuffer<ByteBuffer> args) {
			getCLKernel()
				.setArg(0, energies)
				.setArg(1, numEnergies)
				.setArg(2, efuncArgs)
				.setArg(3, args)
				.setNullArg(4, (EnergyReadThreads + 1)*Double.BYTES);
		}
		
		public void runAsync() {
			super.runAsync(EnergyReadThreads, EnergyReadThreads);
		}
	}
}
