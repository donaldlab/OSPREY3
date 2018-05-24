/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.List;

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
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
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
		
		public DofInfo(FreeDihedral dof, BigForcefieldEnergy.Subset subset) {
			this.res = dof.getResidue();
			this.subset = subset;
			this.dihedralIndices = res.template.getDihedralDefiningAtoms(dof.getDihedralNumber());
			this.rotatedIndices = res.template.getDihedralRotatedAtoms(dof.getDihedralNumber());
			this.atomOffset = ffenergy.getAtomOffset(res);
			this.numAtoms = res.atoms.size();
		}
	}
	
	private BigForcefieldEnergy ffenergy;
	private int ffSequenceNumber;
	
	private Kernel.Function func;
	private static Integer blockThreads = null;
	
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
	
	public CCDKernelCuda(GpuStream stream)
	throws IOException {
		super(stream, "ccd");
		
		ffargs = stream.makeByteBuffer(48);
	}
	
	@Deprecated
	public void init(MoleculeModifierAndScorer mof) {
		init(new MoleculeObjectiveFunction(mof));
	}
	
	public void init(MoleculeObjectiveFunction mof) {
		
		GpuStream stream = getStream();
		
		// make sure this thread can use the cuda context
		stream.getContext().attachCurrentThread();
		
		// get the energy function
		if (mof.efunc instanceof BigForcefieldEnergy) {
			this.ffenergy = (BigForcefieldEnergy)mof.efunc;
		} else {
			throw new Error("CCD kernel needs a " + BigForcefieldEnergy.class.getSimpleName() + ", not a " + mof.efunc.getClass().getSimpleName() + ". this is a bug.");
		}
		
		// handle any chemical changes
		ffSequenceNumber = ffenergy.getFullSubset().handleChemicalChanges();
		ffenergy.updateCoords();
		
		// wrap the incoming buffers
		coords = stream.makeOrExpandBuffer(coords, ffenergy.getCoords());
		atomFlags = stream.makeOrExpandBuffer(atomFlags, ffenergy.getAtomFlags());
		precomputed = stream.makeOrExpandBuffer(precomputed, ffenergy.getPrecomputed());
		
		// make the args buffer
		ByteBuffer argsBuf = ffargs.getHostBuffer();
		argsBuf.rewind();
		argsBuf.putInt(ffenergy.getFullSubset().getNumAtomPairs());
		argsBuf.putInt(ffenergy.getFullSubset().getNum14AtomPairs());
		argsBuf.putDouble(ffenergy.getParams().coulombFactor);
		argsBuf.putDouble(ffenergy.getParams().scaledCoulombFactor);
		argsBuf.putDouble(ffenergy.getParams().solvationCutoff2);
		argsBuf.putDouble(ffenergy.getFullSubset().getInternalSolvationEnergy());
		argsBuf.put((byte)(ffenergy.getParams().useDistDependentDielectric ? 1 : 0));
		argsBuf.put((byte)(ffenergy.getParams().useHElectrostatics ? 1 : 0));
		argsBuf.put((byte)(ffenergy.getParams().useHVdw ? 1 : 0));
		argsBuf.put((byte)(ffenergy.getParams().useEEF1 ? 1 : 0));
		
		// upload static forcefield info
		atomFlags.uploadAsync();
		precomputed.uploadAsync();
		ffargs.uploadAsync();
		
		// get info about the dofs
		dofInfos = new ArrayList<>();
		int subsetsSize = 0;
		int numRotatedAtoms = 0;
		for (int d=0; d<mof.getNumDOFs(); d++) {
			
			// make sure the DoF is a FreeDihedral. that's all we support on the gpu at the moment
			DegreeOfFreedom dofBase = mof.pmol.dofs.get(d);
			if (!(dofBase instanceof FreeDihedral)) {
				throw new Error("degree-of-freedom type " + dofBase.getClass().getSimpleName() + " not yet supported by CCD kernel."
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
		
		// allocate dof-related buffers
		subsetTables = stream.makeOrExpandIntBuffer(subsetTables, subsetsSize);
		IntBuffer subsetTablesBuf = subsetTables.getHostBuffer();
		subsetTablesBuf.clear();
		
		dofargs = stream.makeOrExpandByteBuffer(dofargs, dofInfos.size()*48);
		ByteBuffer dofargsBuf = dofargs.getHostBuffer();
		dofargsBuf.clear();
		
		dihedralIndices = stream.makeOrExpandIntBuffer(dihedralIndices, dofInfos.size()*4);
		IntBuffer dihedralIndicesBuf = dihedralIndices.getHostBuffer();
		dihedralIndicesBuf.clear();
		
		rotatedIndices = stream.makeOrExpandIntBuffer(rotatedIndices, numRotatedAtoms);
		IntBuffer rotatedIndicesBuf = rotatedIndices.getHostBuffer();
		rotatedIndicesBuf.clear();
		
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
			dofargsBuf.putDouble(dofInfo.subset.getInternalSolvationEnergy());
			
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
		final int fMaxNumCoords = maxNumCoords;
		
		// upload more bufs
		dofargs.uploadAsync();
		subsetTables.uploadAsync();
		dihedralIndices.uploadAsync();
		rotatedIndices.uploadAsync();
		
		// init other ccd inputs,outputs
		xAndBounds = stream.makeOrExpandDoubleBuffer(xAndBounds, dofInfos.size()*3);
		ccdOut = stream.makeOrExpandDoubleBuffer(ccdOut, dofInfos.size() + 1);
		
		// init the kernel function
		func = makeFunction("ccd");
		func.numBlocks = 1;
		func.sharedMemCalc = new Kernel.SharedMemCalculator() {
			@Override
			public int calcBytes(int blockThreads) {
				return blockThreads*Double.BYTES // energy reduction
					+ fMaxNumCoords*Double.BYTES // coords copy
					+ dofInfos.size()*Double.BYTES // nextx
					+ dofInfos.size()*Double.BYTES // firstSteps
					+ dofInfos.size()*Double.BYTES; // lastSteps
			}
		};
		func.setArgs(Pointer.to(
			coords.getDevicePointer(),
			atomFlags.getDevicePointer(),
			precomputed.getDevicePointer(),
			ffargs.getDevicePointer(),
			subsetTables.getDevicePointer(),
			dihedralIndices.getDevicePointer(),
			rotatedIndices.getDevicePointer(),
			dofargs.getDevicePointer(),
			Pointer.to(new int[] { maxNumCoords }),
			xAndBounds.getDevicePointer(),
			Pointer.to(new int[] { dofInfos.size() }),
			ccdOut.getDevicePointer()
		));
		
		// calc the number of block threads
		if (blockThreads == null) {
			blockThreads = func.calcMaxBlockThreads();
		}
		func.blockThreads = blockThreads;
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
		result.energy = buf.get();
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
		if (coords != null) {
			
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
		
			coords = null;
		}
	}
	
	@Override
	protected void finalize()
	throws Throwable {
		try {
			if (coords != null) {
				System.err.println("WARNING: " + getClass().getName() + " was garbage collected, but not cleaned up. Attempting cleanup now");
				cleanup();
			}
		} finally {
			super.finalize();
		}
	}
		
}
