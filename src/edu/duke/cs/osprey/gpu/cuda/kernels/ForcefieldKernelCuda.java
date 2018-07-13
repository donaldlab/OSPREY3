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

import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.gpu.ForcefieldKernel;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.tools.MathTools;
import jcuda.Pointer;

public class ForcefieldKernelCuda extends Kernel implements ForcefieldKernel {
	
	private Kernel.Function func;
	
	private CUBuffer<DoubleBuffer> coords;
	private CUBuffer<IntBuffer> atomFlags;
	private CUBuffer<DoubleBuffer> precomputed;
	private CUBuffer<IntBuffer> subsetTable;
	private CUBuffer<DoubleBuffer> energies;
	
	private CUBuffer<ByteBuffer> args;
	
	private BigForcefieldEnergy ffenergy;
	private BigForcefieldEnergy.Subset subset;
	
	public ForcefieldKernelCuda(GpuStream stream, BigForcefieldEnergy ffenergy)
	throws IOException {
		super(stream, "forcefield");
		
		this.ffenergy = ffenergy;
		
		func = makeFunction("calc");
		func.blockThreads = 512; // NOTE: optimizing this doesn't make much difference
		func.sharedMemCalc = new Kernel.SharedMemCalculator() {
			@Override
			public int calcBytes(int blockThreads) {
				return blockThreads*Double.BYTES;
			}
		};
		
		// allocate the buffers
		coords = getStream().makeBuffer(ffenergy.getCoords());
		atomFlags = getStream().makeBuffer(ffenergy.getAtomFlags());
		precomputed = getStream().makeBuffer(ffenergy.getPrecomputed());
		subsetTable = getStream().makeIntBuffer(ffenergy.getFullSubset().getNumAtomPairs());
		energies = getStream().makeDoubleBuffer(getEnergySize(ffenergy.getFullSubset(), func.blockThreads));
		
		// upload static info
		atomFlags.uploadAsync();
		precomputed.uploadAsync();
		
		// make the args buffer
		args = getStream().makeByteBuffer(40);
		ByteBuffer argsBuf = args.getHostBuffer();
		argsBuf.rewind();
		argsBuf.putInt(0); // set by setSubsetInternal()
		argsBuf.putInt(0); // 
		argsBuf.putDouble(ffenergy.getParams().coulombFactor);
		argsBuf.putDouble(ffenergy.getParams().scaledCoulombFactor);
		argsBuf.putDouble(ffenergy.getParams().solvationCutoff2);
		argsBuf.put((byte)(ffenergy.getParams().useDistDependentDielectric ? 1 : 0));
		argsBuf.put((byte)(ffenergy.getParams().useHElectrostatics ? 1 : 0));
		argsBuf.put((byte)(ffenergy.getParams().useHVdw ? 1 : 0));
		argsBuf.put((byte)0); // set by setSubsetInternal()
		argsBuf.put((byte)(ffenergy.getParams().useEEF1 ? 1 : 0));
		argsBuf.flip();
		
		// set the subset
		// NOTE: setting the subset uploads the args too
		subset = null;
		setSubsetInternal(ffenergy.getFullSubset());
		
		func.setArgs(Pointer.to(
			coords.getDevicePointer(),
			atomFlags.getDevicePointer(),
			precomputed.getDevicePointer(),
			subsetTable.getDevicePointer(),
			args.getDevicePointer(),
			energies.getDevicePointer()
		));
	}
	
	public CUBuffer<DoubleBuffer> getCoords() {
		return coords;
	}
	
	public CUBuffer<IntBuffer> getAtomFlags() {
		return atomFlags;
	}
	
	public CUBuffer<DoubleBuffer> getPrecomputed() {
		return precomputed;
	}
	
	public CUBuffer<IntBuffer> getSubsetTable() {
		return subsetTable;
	}
	
	public CUBuffer<DoubleBuffer> getEnergies() {
		return energies;
	}
	
	public CUBuffer<ByteBuffer> getArgs() {
		return args;
	}
	
	@Override
	public BigForcefieldEnergy getForcefield() {
		return ffenergy;
	}
	
	@Override
	public BigForcefieldEnergy.Subset getSubset() {
		return subset;
	}
	
	@Override
	public boolean setSubset(BigForcefieldEnergy.Subset subset) {
		return setSubsetInternal(subset);
	}
	
	private boolean setSubsetInternal(BigForcefieldEnergy.Subset subset) {
		
		// short circuit: don't change things unless we need to
		if (this.subset == subset) {
			return false;
		}
		
		this.subset = subset;
		boolean useSubset = subset.getSubsetTable() != null;
		
		func.numBlocks = MathTools.divUp(subset.getNumAtomPairs(), func.blockThreads);
		
		// make sure this thread can use the cuda context
		getStream().getContext().attachCurrentThread();
		
		// update kernel args and upload
		ByteBuffer buf = args.getHostBuffer();
		buf.putInt(0, subset.getNumAtomPairs());
		buf.putInt(4, subset.getNum14AtomPairs());
		buf.put(35, (byte)(useSubset ? 1 : 0));
		buf.rewind();
		args.uploadAsync();
		
		if (useSubset) {
			
			// upload subset table
			subsetTable.getHostBuffer().clear();
			subset.getSubsetTable().rewind();
			subsetTable.getHostBuffer().put(subset.getSubsetTable());
			subsetTable.getHostBuffer().flip();
			subsetTable.uploadAsync();
		}
		
		return true;
	}
	
	@Override
	public void runAsync() {
		
		// make sure this thread can use the cuda context
		getStream().getContext().attachCurrentThread();
		
		func.runAsync();
	}
	
	private static int getEnergySize(BigForcefieldEnergy.Subset subset, int blockThreads) {
		return MathTools.divUp(subset.getNumAtomPairs(), blockThreads);
	}
	
	@Override
	public void cleanup() {
		if (coords != null) {
			
			coords.cleanup();
			atomFlags.cleanup();
			precomputed.cleanup();
			subsetTable.cleanup();
			args.cleanup();
			energies.cleanup();
			
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
	
	@Override
	public void uploadCoordsAsync() {
		
		// make sure this thread can use the cuda context
		getStream().getContext().attachCurrentThread();
		
		// tell the forcefield to gather updated coords
		ffenergy.updateCoords();
		
		coords.uploadAsync();
	}

	@Override
	public double downloadEnergySync() {
		
		// make sure this thread can use the cuda context
		getStream().getContext().attachCurrentThread();
		
		energies.downloadSync();
		DoubleBuffer buf = energies.getHostBuffer();
		
		// do the last bit of the energy sum on the cpu
		// add one element per work group on the gpu
		// typically, it's a factor of groupSize less than the number of atom pairs
		double energy = subset.getInternalSolvationEnergy();
		buf.rewind();
		int n = getEnergySize(subset, func.blockThreads);
		for (int i=0; i<n; i++) {
			energy += buf.get();
		}
		return energy;
	}
}
