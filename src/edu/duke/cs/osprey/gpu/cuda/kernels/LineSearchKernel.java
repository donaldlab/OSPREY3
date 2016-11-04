package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.List;

import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.structure.Residue;
import jcuda.Pointer;
import jcuda.driver.CUlimit;
import jcuda.driver.JCudaDriver;

public class LineSearchKernel extends Kernel {
	
	private CUBuffer<IntBuffer> dihedralIndices;
	private CUBuffer<IntBuffer> rotatedIndices;
	private CUBuffer<DoubleBuffer> lsargs;
	
	private CUBuffer<DoubleBuffer> result;
	
	private Pointer pKernelArgs;
	
	private ForcefieldKernelCuda ffKernel;
	
	public LineSearchKernel(ForcefieldKernelCuda ffKernel, FreeDihedral dof)
	throws IOException {
		super(ffKernel.getStream(), "linesearch", "calc");
		
		this.ffKernel = ffKernel;
		
		// get the dof indices
		Residue res = dof.getResidue();
		int coordsOffset = ffKernel.getForcefield().getAtomOffset(res);
	
		int[] dihedralIndicesSrc = res.template.getDihedralDefiningAtoms(dof.getDihedralNumber());
		dihedralIndices = getStream().makeIntBuffer(dihedralIndicesSrc.length);
		for (int i=0; i<dihedralIndicesSrc.length; i++) {
			dihedralIndices.getHostBuffer().put(dihedralIndicesSrc[i] + coordsOffset);
		}
		dihedralIndices.getHostBuffer().flip();
		dihedralIndices.uploadAsync();
		
		List<Integer> rotatedIndicesSrc = res.template.getDihedralRotatedAtoms(dof.getDihedralNumber());
		rotatedIndices = getStream().makeIntBuffer(rotatedIndicesSrc.size());
		for (int i=0; i<rotatedIndicesSrc.size(); i++) {
			rotatedIndices.getHostBuffer().put(rotatedIndicesSrc.get(i) + coordsOffset);
		}
		rotatedIndices.getHostBuffer().flip();
		rotatedIndices.uploadAsync();
		
		// allocate the rest of the buffers
		lsargs = getStream().makeDoubleBuffer(4);
		result = getStream().makeDoubleBuffer(2);
		
		pKernelArgs = Pointer.to(
			ffKernel.getCoords().makeDevicePointer(),
			dihedralIndices.makeDevicePointer(),
			Pointer.to(new int[] { rotatedIndicesSrc.size() }),
			rotatedIndices.makeDevicePointer(),
			ffKernel.getAtomFlags().makeDevicePointer(),
			ffKernel.getPrecomputed().makeDevicePointer(),
			ffKernel.getSubsetTable().makeDevicePointer(),
			ffKernel.getArgs().makeDevicePointer(),
			lsargs.makeDevicePointer(),
			result.makeDevicePointer()
		);
		
		// TEMP: try optimizing dynamic parallelism options
		JCudaDriver.cuCtxSetLimit(CUlimit.CU_LIMIT_DEV_RUNTIME_SYNC_DEPTH, 1);
		JCudaDriver.cuCtxSetLimit(CUlimit.CU_LIMIT_DEV_RUNTIME_PENDING_LAUNCH_COUNT, 1);
	}
	
	public void runAsync(double xdmin, double xdmax, double xd, double step) {
		
		// upload args
		DoubleBuffer lsargsBuf = lsargs.getHostBuffer();
		lsargsBuf.rewind();
		lsargsBuf.put(xdmin);
		lsargsBuf.put(xdmax);
		lsargsBuf.put(xd);
		lsargsBuf.put(step);
		lsargs.uploadAsync();
		
		runAsync(1, 1, 0, pKernelArgs);
	}
	
	public void cleanup() {
		dihedralIndices.cleanup();
		rotatedIndices.cleanup();
		result.cleanup();
	}
	
	public void downloadResultSync() {
		result.downloadSync();
	}
	
	public double getXdstar() {
		return result.getHostBuffer().get(0);
	}
	
	public double getFxdstar() {
		return result.getHostBuffer().get(1) + ffKernel.getSubset().getInternalSolvationEnergy();
	}
}
