package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.structure.Residue;
import jcuda.Pointer;

public class CCDKernelCuda extends Kernel {
	
	private class DofInfo {
		
		public final Residue res;
		public final BigForcefieldEnergy.Subset subset;
		public final int[] dihedralIndices;
		public final List<Integer> rotatedIndices;
		public final int atomOffset;
		public final int numAtoms;
		
		public DofInfo(FreeDihedral dof) {
			
			this.res = dof.getResidue();
			this.subset = ffenergy.new Subset(ffenergy.getInteractions().makeSubsetByResidue(res));
			
			this.dihedralIndices = res.template.getDihedralDefiningAtoms(dof.getDihedralNumber());
			this.rotatedIndices = res.template.getDihedralRotatedAtoms(dof.getDihedralNumber());
			this.atomOffset = ffenergy.getAtomOffset(res);
			this.numAtoms = res.atoms.size();
		}
	}
	
	public static class Result {
		
		public final DoubleMatrix1D x;
		public final double energy;
		
		public Result(int numDofs, DoubleBuffer buf, double internalSolvationEnergy) {
			x = DoubleFactory1D.dense.make(numDofs);
			buf.rewind();
			energy = buf.get() + internalSolvationEnergy;
			for (int d=0; d<numDofs; d++) {
				x.set(d, buf.get());
			}
		}
	}
	
	// OPTIMIZATION: generally, more threads = more faster, but more threads use more GPU SM resources
	// this default value will probably under-saturate newer cards and require too many resources for older cards
	// ideally, this should be optimized for each hardware platform
	private static final int DefaultNumThreads = 512;
		
	
	private BigForcefieldEnergy ffenergy;
	
	private Kernel.Function func;
	
	private CUBuffer<DoubleBuffer> coords;
	private CUBuffer<IntBuffer> atomFlags;
	private CUBuffer<DoubleBuffer> precomputed;
	private CUBuffer<ByteBuffer> ffargs;
	
	private CUBuffer<ByteBuffer> dofargs;
	private CUBuffer<IntBuffer> subsetTables;
	private CUBuffer<IntBuffer> dihedralIndices;
	private CUBuffer<IntBuffer> rotatedIndices;
	
	private CUBuffer<DoubleBuffer> xAndBounds;
	private CUBuffer<DoubleBuffer> ccdOut;
	
	private List<DofInfo> dofInfos;
	
	public CCDKernelCuda(GpuStream stream, BigForcefieldEnergy ffenergy, List<FreeDihedral> dofs)
	throws IOException {
		this(stream, ffenergy, dofs, DefaultNumThreads);
	}
	
	public CCDKernelCuda(GpuStream stream, BigForcefieldEnergy ffenergy, List<FreeDihedral> dofs, int numThreads)
	throws IOException {
		super(stream, "ccd");
		
		this.ffenergy = ffenergy;
		
		// wrap the incoming buffers
		coords = stream.makeBuffer(ffenergy.getCoords());
		atomFlags = stream.makeBuffer(ffenergy.getAtomFlags());
		precomputed = stream.makeBuffer(ffenergy.getPrecomputed());
		
		// make the args buffer
		ffargs = stream.makeByteBuffer(40);
		ByteBuffer argsBuf = ffargs.getHostBuffer();
		argsBuf.rewind();
		argsBuf.putInt(ffenergy.getFullSubset().getNumAtomPairs());
		argsBuf.putInt(ffenergy.getFullSubset().getNum14AtomPairs());
		argsBuf.putDouble(ffenergy.getCoulombFactor());
		argsBuf.putDouble(ffenergy.getScaledCoulombFactor());
		argsBuf.putDouble(ffenergy.getSolvationCutoff2());
		argsBuf.put((byte)(ffenergy.useDistDependentDielectric() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHElectrostatics() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHVdw() ? 1 : 0));
		
		// upload static forcefield info
		atomFlags.uploadAsync();
		precomputed.uploadAsync();
		ffargs.uploadAsync();
		
		// get info about the dofs
		dofInfos = new ArrayList<>();
		int subsetsSize = 0;
		int numRotatedAtoms = 0;
		for (FreeDihedral dof : dofs) {
			
			DofInfo dofInfo = new DofInfo(dof);
			dofInfos.add(dofInfo);
			
			// update counts
			subsetsSize += dofInfo.subset.getNumAtomPairs();
			numRotatedAtoms += dofInfo.rotatedIndices.size();
		}
		
		// allocate dof-related buffers
		subsetTables = stream.makeIntBuffer(subsetsSize);
		IntBuffer subsetTablesBuf = subsetTables.getHostBuffer();
		subsetTablesBuf.clear();
		
		dofargs = stream.makeByteBuffer(dofs.size()*40);
		ByteBuffer dofargsBuf = dofargs.getHostBuffer();
		dofargsBuf.clear();
		
		dihedralIndices = stream.makeIntBuffer(dofs.size()*4);
		IntBuffer dihedralIndicesBuf = dihedralIndices.getHostBuffer();
		dihedralIndicesBuf.clear();
		
		rotatedIndices = stream.makeIntBuffer(numRotatedAtoms);
		IntBuffer rotatedIndicesBuf = rotatedIndices.getHostBuffer();
		rotatedIndicesBuf.clear();
		
		// populate dof buffers
		int maxNumCoords = 0;
		for (int d=0; d<dofs.size(); d++) {
			DofInfo dofInfo = dofInfos.get(d);
			
			int firstCoord = dofInfo.atomOffset*3;
			int lastCoord = (dofInfo.atomOffset + dofInfo.numAtoms)*3 - 1;
			int numCoords = lastCoord - firstCoord + 1;
			maxNumCoords = Math.max(maxNumCoords, numCoords);
			
			// update dofargs
			dofargsBuf.putInt(subsetTablesBuf.position());
			dofargsBuf.putInt(dofInfo.subset.getNumAtomPairs());
			dofargsBuf.putInt(dofInfo.subset.getNum14AtomPairs());
			dofargsBuf.putInt(rotatedIndicesBuf.position());
			dofargsBuf.putInt(dofInfo.rotatedIndices.size());
			dofargsBuf.putInt(firstCoord);
			dofargsBuf.putInt(lastCoord);
			dofargsBuf.putInt(0); // spacer for word alignment
			
			// append subset table
			dofInfo.subset.getSubsetTable().rewind();
			subsetTablesBuf.put(dofInfo.subset.getSubsetTable());
			
			// collect dihedral indices
			for (int i=0; i<dofInfo.dihedralIndices.length; i++) {
				dihedralIndicesBuf.put(dofInfo.dihedralIndices[i]);
			}
			
			// collect rotated indices
			for (int i=0; i<dofInfo.rotatedIndices.size(); i++) {
				rotatedIndicesBuf.put(dofInfo.rotatedIndices.get(i));
			}
		}
		
		// upload more bufs
		dofargs.uploadAsync();
		subsetTables.uploadAsync();
		dihedralIndices.uploadAsync();
		rotatedIndices.uploadAsync();
		
		// init other ccd inputs,outputs
		xAndBounds = stream.makeDoubleBuffer(dofs.size()*3);
		ccdOut = stream.makeDoubleBuffer(dofs.size() + 1);
		
		// init the kernel function
		func = makeFunction("ccd");
		func.numBlocks = 1;
		func.blockThreads = numThreads;
		func.sharedMemBytes = numThreads*Double.BYTES // energy reduction
			+ maxNumCoords*Double.BYTES // coords copy
			+ dofs.size()*Double.BYTES // nextx
			+ dofs.size()*Double.BYTES // firstSteps
			+ dofs.size()*Double.BYTES; // lastSteps
		func.setArgs(Pointer.to(
			coords.makeDevicePointer(),
			atomFlags.makeDevicePointer(),
			precomputed.makeDevicePointer(),
			ffargs.makeDevicePointer(),
			subsetTables.makeDevicePointer(),
			dihedralIndices.makeDevicePointer(),
			rotatedIndices.makeDevicePointer(),
			dofargs.makeDevicePointer(),
			Pointer.to(new int[] { maxNumCoords }),
			xAndBounds.makeDevicePointer(),
			Pointer.to(new int[] { dofs.size() }),
			ccdOut.makeDevicePointer()
		));
	}
	
	public void uploadCoordsAsync() {
		
		// tell the forcefield to gather updated coords
		ffenergy.updateCoords();
		
		coords.uploadAsync();
	}
	
	public void runAsync(DoubleMatrix1D x, ObjectiveFunction.DofBounds dofBounds) {
		
		// upload x and bounds
		DoubleBuffer buf = xAndBounds.getHostBuffer();
		buf.clear();
		int numDofs = dofInfos.size();
		for (int d=0; d<numDofs; d++) {
			buf.put(Math.toRadians(x.get(d)));
			buf.put(Math.toRadians(dofBounds.getMin(d)));
			buf.put(Math.toRadians(dofBounds.getMax(d)));
		}
		xAndBounds.uploadAsync();
		
		// run the kernel
		func.runAsync();
	}
	
	public Result downloadResultSync() {
		return new Result(
			dofInfos.size(),
			ccdOut.downloadSync(),
			ffenergy.getFullSubset().getInternalSolvationEnergy()
		);
	}
	
	public Result runSync(DoubleMatrix1D x, ObjectiveFunction.DofBounds dofBounds) {
		runAsync(x, dofBounds);
		return downloadResultSync();
	}
	
	public void cleanup() {
		coords.cleanup();
		atomFlags.cleanup();
		precomputed.cleanup();
		ffargs.cleanup();
		dofargs.cleanup();
		subsetTables.cleanup();
		dihedralIndices.cleanup();
		rotatedIndices.cleanup();
		xAndBounds.cleanup();
		ccdOut.cleanup();
	}
}
