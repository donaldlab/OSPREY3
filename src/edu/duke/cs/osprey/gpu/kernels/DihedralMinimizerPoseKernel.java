package edu.duke.cs.osprey.gpu.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLMemory;

import edu.duke.cs.osprey.gpu.Kernel;
import edu.duke.cs.osprey.gpu.GpuQueue;

public class DihedralMinimizerPoseKernel extends Kernel {
	
	private int numRotatedAtoms;
	private CLBuffer<IntBuffer> dihedralIndices;
	private CLBuffer<IntBuffer> rotatedIndices;
	private CLBuffer<ByteBuffer> args;
	
	public DihedralMinimizerPoseKernel(GpuQueue queue)
	throws IOException {
		super(queue, "dihedralMinimizer.cl", "pose");
	}
	
	public void init(CLBuffer<DoubleBuffer> coords, int[] dihedralIndices, int[] rotatedIndices) {
		
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
		args = context.createByteBuffer(16, CLMemory.Mem.READ_WRITE);
		ByteBuffer buf = args.getBuffer();
		buf.rewind();
		buf.putInt(rotatedIndices.length);
		buf.putInt(0); // 4-byte pad for 8-byte alignment
		buf.putDouble(0); // set later by setDihedral, TODO: use initial dihedral?
		buf.flip();
		uploadBufferAsync(args);
		
		// set kernel args
		getCLKernel()
			.setArg(0, coords)
			.setArg(1, this.dihedralIndices)
			.setArg(2, this.rotatedIndices)
			.setArg(3, args);
		
		numRotatedAtoms = rotatedIndices.length;
	}
	
	public void setDihedralAsync(double val) {
		ByteBuffer buf = args.getBuffer();
		buf.putDouble(8, val);
		buf.rewind();
		uploadBufferAsync(args);
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
