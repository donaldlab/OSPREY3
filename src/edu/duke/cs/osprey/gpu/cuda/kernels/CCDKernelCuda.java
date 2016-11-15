package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Profiler;
import jcuda.Pointer;

public class CCDKernelCuda extends Kernel {
	
	private class DofInfo {
		
		public final Residue res;
		public final BigForcefieldEnergy.Subset subset;
		public final int[] dihedralIndices;
		public final List<Integer> rotatedIndices;
		public final int atomOffset;
		public final int numAtoms;
		
		public DofInfo(FreeDihedral dof, BigForcefieldEnergy.Subset subset) {
			this.res = dof.getResidue();
			this.subset = subset;
			this.dihedralIndices = res.template.getDihedralDefiningAtoms(dof.getDihedralNumber());
			this.rotatedIndices = res.template.getDihedralRotatedAtoms(dof.getDihedralNumber());
			this.atomOffset = ffenergy.getAtomOffset(res);
			this.numAtoms = res.atoms.size();
		}
	}
	
	// OPTIMIZATION: generally, more threads = more faster, but more threads use more GPU SM resources
	// this default value will probably under-saturate newer cards and require too many resources for older cards
	// ideally, this should be optimized for each hardware platform
	private static final int DefaultNumThreads = 512;
	
	private BigForcefieldEnergy ffenergy;
	private int ffSequenceNumber;
	
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
	
	public CCDKernelCuda(GpuStream stream, MoleculeModifierAndScorer mof)
	throws IOException {
		this(stream, mof, DefaultNumThreads);
	}
	
	public CCDKernelCuda(GpuStream stream, MoleculeModifierAndScorer mof, int numThreads)
	throws IOException {
		super(stream, "ccd");
		
		// TEMP
		Profiler profiler = new Profiler();
		profiler.start("efunc");
		
		// get the energy function
		if (mof.getEfunc() instanceof BigForcefieldEnergy) {
			this.ffenergy = (BigForcefieldEnergy)mof.getEfunc();
		} else {
			throw new Error("CCD kernel needs a " + BigForcefieldEnergy.class.getSimpleName() + ", not a " + mof.getEfunc().getClass().getSimpleName() + ". this is a bug.");
		}
		ffSequenceNumber = ffenergy.getFullSubset().handleChemicalChanges();
		
		// TEMP
		profiler.start("ff allocate"); // SLOW
		
		// wrap the incoming buffers
		coords = stream.makeBuffer(ffenergy.getCoords());
		atomFlags = stream.makeBuffer(ffenergy.getAtomFlags());
		precomputed = stream.makeBuffer(ffenergy.getPrecomputed());
		
		// TEMP
		profiler.start("args");
		
		// make the args buffer
		ffargs = stream.makeByteBuffer(40);
		ByteBuffer argsBuf = ffargs.getHostBuffer();
		argsBuf.rewind();
		argsBuf.putInt(ffenergy.getFullSubset().getNumAtomPairs());
		argsBuf.putInt(ffenergy.getFullSubset().getNum14AtomPairs());
		argsBuf.putDouble(ffenergy.getParams().coulombFactor);
		argsBuf.putDouble(ffenergy.getParams().scaledCoulombFactor);
		argsBuf.putDouble(ffenergy.getParams().solvationCutoff2);
		argsBuf.put((byte)(ffenergy.getParams().useDistDependentDielectric ? 1 : 0));
		argsBuf.put((byte)(ffenergy.getParams().useHElectrostatics ? 1 : 0));
		argsBuf.put((byte)(ffenergy.getParams().useHVdw ? 1 : 0));
		
		// TEMP
		profiler.start("ff upload");
		
		// upload static forcefield info
		atomFlags.uploadAsync();
		precomputed.uploadAsync();
		ffargs.uploadAsync();
		
		// TEMP
		profiler.start("dofs");
		
		// get info about the dofs
		dofInfos = new ArrayList<>();
		int subsetsSize = 0;
		int numRotatedAtoms = 0;
		for (int d=0; d<mof.getNumDOFs(); d++) {
			
			// make sure the DoF is a FreeDihedral. that's all we support on the gpu at the moment
			DegreeOfFreedom dofBase = mof.getDOFs().get(d);
			if (!(dofBase instanceof FreeDihedral)) {
				throw new Error("degree-of-freedom type " + dofBase.getClass().getSimpleName() + " not yet supported by CCD kerne."
					+ " Use CPU minimizer with GPU energy function instead");
			}
			
			// get the dof and its ff subset
			FreeDihedral dof = (FreeDihedral)dofBase;
			BigForcefieldEnergy.Subset subset = (BigForcefieldEnergy.Subset)mof.getEfunc(d);
			
			DofInfo dofInfo = new DofInfo(dof, subset);
			dofInfos.add(dofInfo);
			
			// update counts
			subsetsSize += dofInfo.subset.getNumAtomPairs();
			numRotatedAtoms += dofInfo.rotatedIndices.size();
		}
		
		// TEMP
		profiler.start("dof allocate");
		
		// allocate dof-related buffers
		subsetTables = stream.makeIntBuffer(subsetsSize);
		IntBuffer subsetTablesBuf = subsetTables.getHostBuffer();
		subsetTablesBuf.clear();
		
		dofargs = stream.makeByteBuffer(dofInfos.size()*40);
		ByteBuffer dofargsBuf = dofargs.getHostBuffer();
		dofargsBuf.clear();
		
		dihedralIndices = stream.makeIntBuffer(dofInfos.size()*4);
		IntBuffer dihedralIndicesBuf = dihedralIndices.getHostBuffer();
		dihedralIndicesBuf.clear();
		
		rotatedIndices = stream.makeIntBuffer(numRotatedAtoms);
		IntBuffer rotatedIndicesBuf = rotatedIndices.getHostBuffer();
		rotatedIndicesBuf.clear();
		
		// TEMP
		profiler.start("dof populate");
		
		// populate dof buffers
		int maxNumCoords = 0;
		for (int d=0; d<dofInfos.size(); d++) {
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
		
		// TEMP
		profiler.start("dof upload");
		
		// upload more bufs
		dofargs.uploadAsync();
		subsetTables.uploadAsync();
		dihedralIndices.uploadAsync();
		rotatedIndices.uploadAsync();
		
		// TEMP
		profiler.start("func");
		
		// init other ccd inputs,outputs
		xAndBounds = stream.makeDoubleBuffer(dofInfos.size()*3);
		ccdOut = stream.makeDoubleBuffer(dofInfos.size() + 1);
		
		// init the kernel function
		func = makeFunction("ccd");
		func.numBlocks = 1;
		func.blockThreads = numThreads;
		func.sharedMemBytes = numThreads*Double.BYTES // energy reduction
			+ maxNumCoords*Double.BYTES // coords copy
			+ dofInfos.size()*Double.BYTES // nextx
			+ dofInfos.size()*Double.BYTES // firstSteps
			+ dofInfos.size()*Double.BYTES; // lastSteps
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
			Pointer.to(new int[] { dofInfos.size() }),
			ccdOut.makeDevicePointer()
		));
		
		// TEMP
		profiler.stop();
		System.out.println("\tccd kernel: " + profiler.makeReport(TimeUnit.MILLISECONDS));
	}
	
	public void uploadCoordsAsync() {
		
		// tell the forcefield to gather updated coords
		ffenergy.updateCoords();
		
		coords.uploadAsync();
	}
	
	public void runAsync(DoubleMatrix1D x, ObjectiveFunction.DofBounds dofBounds) {
		
		// make sure the energy function is still valid
		if (ffenergy.getFullSubset().handleChemicalChanges() != ffSequenceNumber) {
			throw new Error("don't re-use kernel instances after chemical changes. This is a bug");
		}
		
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
	
	public Minimizer.Result downloadResultSync() {
		
		// allocate space for the result
		int numDofs = dofInfos.size();
		Minimizer.Result result = new Minimizer.Result(DoubleFactory1D.dense.make(numDofs), 0);
		
		// wait for the kernel to finish and download the out buffer
		DoubleBuffer buf = ccdOut.downloadSync();
		
		// copy to the result
		buf.rewind();
		result.energy = buf.get() + ffenergy.getFullSubset().getInternalSolvationEnergy();
		for (int d=0; d<numDofs; d++) {
			result.dofValues.set(d, Math.toDegrees(buf.get()));
		}
		
		return result;
	}
	
	public Minimizer.Result runSync(DoubleMatrix1D x, ObjectiveFunction.DofBounds dofBounds) {
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
