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

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache.ResPair;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import jcuda.Pointer;

public class ResidueForcefieldEnergyCuda extends Kernel implements EnergyFunction.DecomposableByDof, EnergyFunction.NeedsCleanup {
	
	private static final long serialVersionUID = 4015880661919715967L;
	
	public final GpuStreamPool streams;
	public final ForcefieldParams ffparams;
	public final ResidueInteractions inters;
	public final Residues residues;
	
	public final boolean isBroken;
	
	/* buffer layout:
	 * NOTE: try to use 8-byte alignments to be happy on 64-bit machines
	 * 
	 * long flags (useHEs, useHvdW, distDepDielect, useEEF1)
	 * double coulombFactor
	 * double scaledCoulombFactor
	 * 
	 * for each residue pair:
	 *    long offset
	 * 
	 * for each residue pair:
	 *    long numAtomPairs
	 *    long atomsOffset1
	 *    long atomsOffset2
	 *    double weight
	 *    double offset
	 *    
	 *    // NOTE: use struct-of-arrays here so GPU memory accesses are coalesced
	 *    for each atom pair:
	 *       long flags  (bit isHeavyPair, bit is14Bonded, 6 bits space, 3 byte space, short atomOffset1, short atomOffset2)
	 *    for each atom pair:
	 *       double charge
	 *    for each atom pair:
	 *       double Aij
	 *    for each atom pair:
	 *       double Bij
	 *       
	 *    // if EEF1 is used
	 *    for each atom pair:
	 *       double radius1
	 *    for each atom pair:
	 *       double lambda1
	 *    for each atom pair:
	 *       double alpha1
	 *    for each atom pair:
	 *       double radius2
	 *    for each atom pair:
	 *       double lambda2
	 *    for each atom pair:
	 *       double alpha2
	 */
	private CUBuffer<ByteBuffer> data;
	private CUBuffer<DoubleBuffer> coords;
	private CUBuffer<IntBuffer> allIndices;
	private CUBuffer<DoubleBuffer> energy;
	
	private static final int HeaderBytes = Long.BYTES + Double.BYTES*2;
	private static final int ResPairBytes = Long.BYTES*3 + Double.BYTES*2;
	private static final int AtomPairBytes = Long.BYTES + Double.BYTES*3;
	private static final int EEF1Bytes = Double.BYTES*6;
	
	private ResPair[] resPairs;
	private Map<Residue,Subset> subsets;
	
	private Kernel.Function func;
	private static AtomicInteger blockThreads = new AtomicInteger(-1);
	
	public ResidueForcefieldEnergyCuda(GpuStreamPool streams, ResPairCache resPairCache, ResidueInteractions inters, Molecule mol) {
		this(streams, resPairCache, inters, mol.residues);
	}
	
	public ResidueForcefieldEnergyCuda(GpuStreamPool streams, ResPairCache resPairCache, ResidueInteractions inters, Residues residues) {
		super(streams.checkout(), "residueForcefield");
		
		this.streams = streams;
		this.ffparams = resPairCache.ffparams;
		this.inters = inters;
		this.residues = inters.filter(residues);
		
		// is this a broken conformation?
		for (Residue res : this.residues) {
			if (!res.confProblems.isEmpty()) {
				isBroken = true;
				
				// we're done here, no need to analyze broken conformations
				return;
			}
		}
		isBroken = false;
		
		// compute solvation info if needed
		SolvationForcefield.ResiduesInfo solvInfo = null;
		if (ffparams.solvationForcefield != null) {
			solvInfo = ffparams.solvationForcefield.makeInfo(ffparams, this.residues);
		}
		
		// map the residue numbers to residues
		resPairs = new ResPair[inters.size()];
		int index = 0;
		for (ResidueInteractions.Pair pair : inters) {
			resPairs[index++] = resPairCache.get(this.residues, pair, solvInfo);
		}
		
		subsets = null;
		
		// count atoms and offsets for each residue
		int[] atomOffsetsByResIndex = new int[this.residues.size()];
		Arrays.fill(atomOffsetsByResIndex, -1);
		int atomOffset = 0;
		int numAtoms = 0;
		for (int i=0; i<this.residues.size(); i++) {
			Residue res = this.residues.get(i);
			atomOffsetsByResIndex[i] = atomOffset;
			atomOffset += 3*res.atoms.size();
			numAtoms += res.atoms.size();
		}
		
		GpuStream stream = getStream();
		
		// make the coords buffer
		coords = stream.doubleBuffers.checkout(numAtoms*3);
		
		// allocate the data
		int totalNumAtomPairs = 0;
		for (int i=0; i<resPairs.length; i++) {
			totalNumAtomPairs += resPairs[i].info.numAtomPairs;
		}
		int atomPairBytes = AtomPairBytes + (ffparams.solvationForcefield == SolvationForcefield.EEF1 ? EEF1Bytes : 0);
		data = stream.byteBuffers.checkout(HeaderBytes + (Long.BYTES + ResPairBytes)*resPairs.length + atomPairBytes*totalNumAtomPairs);
		ByteBuffer databuf = data.getHostBuffer();
		
		// put the data header
		long flags = ffparams.hElect ? 1 : 0;
		flags <<= 1;
		flags |= ffparams.hVDW ? 1 : 0;
		flags <<= 1;
		flags |= ffparams.distDepDielect ? 1 : 0;
		flags <<= 1;
		flags |= ffparams.solvationForcefield == SolvationForcefield.EEF1 ? 1 : 0;
		databuf.putLong(flags);
		double coulombFactor = ForcefieldParams.coulombConstant/ffparams.dielectric;
		double scaledCoulombFactor = coulombFactor*ffparams.forcefld.coulombScaling;
		databuf.putDouble(coulombFactor);
		databuf.putDouble(scaledCoulombFactor);
		
		// put the res pair offsets
		long offset = HeaderBytes + Long.BYTES*resPairs.length;
		for (int i=0; i<resPairs.length; i++) {
			databuf.putLong(offset);
			offset += ResPairBytes + atomPairBytes*resPairs[i].info.numAtomPairs;
		}
		
		// put the res pairs and atom pairs
		for (int i=0; i<resPairs.length; i++) {
			ResPair resPair = resPairs[i];
			
			// put the res pair
			databuf.putLong(resPair.info.numAtomPairs);
			databuf.putLong(atomOffsetsByResIndex[resPair.resIndex1]);
			databuf.putLong(atomOffsetsByResIndex[resPair.resIndex2]);
			databuf.putDouble(resPair.weight);
			databuf.putDouble(resPair.offset + resPair.solvEnergy);
			
			// put the atom pairs
			// NOTE: use struct-of-arrays here, not array-of-structs
			// so the GPU can coalesce memory accesses
			for (int j=0; j<resPair.info.numAtomPairs; j++) {
				databuf.putLong(resPair.info.flags[j]);
			}
			for (int k=0; k<resPair.info.numPrecomputedPerAtomPair; k++) {
				for (int j=0; j<resPair.info.numAtomPairs; j++) {
					databuf.putDouble(resPair.info.precomputed[j*resPair.info.numPrecomputedPerAtomPair + k]);
				}
			}
		}
		
		databuf.flip();
		data.uploadAsync();
		
		// make the indices
		allIndices = stream.intBuffers.checkout(resPairs.length);
		IntBuffer allIndicesBuf = allIndices.getHostBuffer();
		for (int i=0; i<resPairs.length; i++) {
			allIndicesBuf.put(i);
		}
		allIndicesBuf.flip();
		allIndices.uploadAsync();
		
		// make the energy buffer
		energy = stream.doubleBuffers.checkout(1);
		
		// make the kernel function
		func = makeFunction("calc");
		func.numBlocks = 1;
		func.sharedMemCalc = (int blockThreads) -> blockThreads*Double.BYTES;
		func.setArgs(Pointer.to(
			coords.getDevicePointer(),
			data.getDevicePointer(),
			Pointer.to(new int[] { 0 }), // NOTE: dummy args here
			Pointer.to(new int[] { 0 }),    // since we're just checking launch resources for now
			energy.getDevicePointer()
		));
		func.blockThreads = func.getBestBlockThreads(blockThreads);
	}
	
	@Override
	public void clean() {
		GpuStream stream = getStream();
		if (coords != null) {
			stream.doubleBuffers.release(coords);
			coords = null;
		}
		if (data != null) {
			stream.byteBuffers.release(data);
			data = null;
		}
		if (allIndices != null) {
			stream.intBuffers.release(allIndices);
			allIndices = null;
		}
		if (energy != null) {
			stream.doubleBuffers.release(energy);
			energy = null;
		}
		if (subsets != null) {
			for (Subset subset : subsets.values()) {
				if (subset.indices != null) {
					stream.intBuffers.release(subset.indices);
				}
			}
			subsets = null;
		}
		streams.release(stream);
	}
	
	@Override
	public double getEnergy() {
		return getEnergy(allIndices);
	}
	
	private double getEnergy(CUBuffer<IntBuffer> indices) {
		
		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}
		
		// make sure this thread can use the cuda context
		getStream().getContext().attachCurrentThread();
		
		// capture the current molecule state
		DoubleBuffer coordsbuf = coords.getHostBuffer();
		coordsbuf.clear();
		for (Residue res : residues) {
			coordsbuf.put(res.coords);
		}
		coordsbuf.clear();
		coords.uploadAsync();
		
		// launch kernel
		func.setArgs(Pointer.to(
			coords.getDevicePointer(),
			data.getDevicePointer(),
			Pointer.to(new int[] { indices.getHostBuffer().limit() }),
			indices.getDevicePointer(),
			energy.getDevicePointer()
		));
		func.runAsync();
		
		// download the energy
		DoubleBuffer buf = energy.downloadSync();
		buf.rewind();
		return buf.get();
	}
	
	private class Subset implements EnergyFunction {
		
		private static final long serialVersionUID = -1749739381007657718L;
		
		private CUBuffer<IntBuffer> indices;
		
		public Subset(Residue res) {
			
			// pass 1: count
			int num = 0;
			for (int i=0; i<resPairs.length; i++) {
				ResPair resPair = resPairs[i];
				if (resPair.res1 == res || resPair.res2 == res) {
					num++;
				}
			}
			
			// no interactions for this residue
			if (num <= 0) {
				return;
			}
			
			// pass 2: collect
			indices = getStream().intBuffers.checkout(num);
			IntBuffer indicesbuf = indices.getHostBuffer();
			for (int i=0; i<resPairs.length; i++) {
				ResPair resPair = resPairs[i];
				if (resPair.res1 == res || resPair.res2 == res) {
					indicesbuf.put(i);
				}
			}
			indicesbuf.flip();
			indices.uploadAsync();
		}
		
		@Override
		public double getEnergy() {
			if (indices == null) {
				return Double.NaN;
			}
			return ResidueForcefieldEnergyCuda.this.getEnergy(indices);
		}
	}
	
	@Override
	public List<EnergyFunction> decomposeByDof(Molecule mol, List<DegreeOfFreedom> dofs) {
		
		if (subsets == null) {
			subsets = new HashMap<>();
		}
		
		List<EnergyFunction> efuncs = new ArrayList<>();
		for (DegreeOfFreedom dof : dofs) {
			Residue res = dof.getResidue();
			
			if (res == null) {
				
				// no res, just use the whole efunc
				efuncs.add(this);
				
			} else {
				
				// make a subset energy function
				Subset subset = subsets.get(res);
				if (subset == null) {
					subset = new Subset(res);
					subsets.put(res, subset);
				}
				efuncs.add(subset);
			}
		}
		
		return efuncs;
	}
}
