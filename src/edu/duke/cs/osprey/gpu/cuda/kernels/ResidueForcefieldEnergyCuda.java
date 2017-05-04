package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
import jcuda.Pointer;

public class ResidueForcefieldEnergyCuda extends Kernel implements EnergyFunction.DecomposableByDof, EnergyFunction.NeedsCleanup {
	
	private static final long serialVersionUID = 4015880661919715967L;
	
	public final GpuStreamPool streams;
	public final ForcefieldParams ffparams;
	public final ResidueInteractions inters;
	
	private boolean isBroken;
	
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
	 *    for each atom pair:
	 *       long flags  (bit isHeavyPair, bit is14Bonded, 6 bits space, 3 byte space, short atomOffset1, short atomOffset2)
	 *       double charge
	 *       double Aij
	 *       double Bij
	 *       double radius1
	 *       double lambda1
	 *       double alpha1
	 *       double radius2
	 *       double lambda2
	 *       double alpha2
	 */
	private CUBuffer<ByteBuffer> data;
	private CUBuffer<DoubleBuffer> coords;
	private CUBuffer<IntBuffer> allIndices;
	private CUBuffer<DoubleBuffer> energy;
	
	private static final int HeaderBytes = Long.BYTES + Double.BYTES*2;
	private static final int ResPairBytes = Long.BYTES*3 + Double.BYTES*2;
	private static final int AtomPairBytes = Short.BYTES*2 + Integer.BYTES + Double.BYTES*9;
	
	private Residue[] ress;
	private ResPair[] resPairs;
	private Map<Residue,Subset> subsets;
	
	private Kernel.Function func;
	
	public ResidueForcefieldEnergyCuda(GpuStreamPool streams, ResPairCache resPairCache, ResidueInteractions inters, Molecule mol)
	throws IOException {
		this(streams, resPairCache, inters, mol.residues);
	}
	
	public ResidueForcefieldEnergyCuda(GpuStreamPool streams, ResPairCache resPairCache, ResidueInteractions inters, List<Residue> residues)
	throws IOException {
		super(streams.checkout(), "residueForcefield");
		
		this.streams = streams;
		this.ffparams = resPairCache.ffparams;
		this.inters = inters;
		
		// compute solvation info if needed
		SolvationForcefield.ResiduesInfo solvInfo = null;
		if (ffparams.solvationForcefield != null) {
			solvInfo = ffparams.solvationForcefield.makeInfo(ffparams, residues);
		}
		
		// TODO: test with no solvation!
		
		// map the residue numbers to residues
		ress = new Residue[inters.getResidueNumbers().size()];
		int index = 0;
		for (String resNum : inters.getResidueNumbers()) {
			ress[index++] = residues.get(Residue.findInOrThrow(residues, resNum));
		}
		resPairs = new ResPair[inters.size()];
		index = 0;
		for (ResidueInteractions.Pair pair : inters) {
			resPairs[index++] = resPairCache.get(ress, pair, solvInfo);
		}
		
		subsets = null;
		
		// is this a broken conformation?
		isBroken = false;
		for (ResPair pair : resPairs) {
			if (pair.res1.confProblems.size() + pair.res2.confProblems.size() > 0) {
				isBroken = true;
				
				// we're done here, no need to analyze broken conformations
				return;
			}
		}
		
		// count atoms and offsets for each residue
		int[] atomOffsetsByResIndex = new int[ress.length];
		Arrays.fill(atomOffsetsByResIndex, -1);
		int atomOffset = 0;
		int numAtoms = 0;
		for (int i=0; i<ress.length; i++) {
			Residue res = ress[i];
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
		data = stream.byteBuffers.checkout(HeaderBytes + (Long.BYTES + ResPairBytes)*resPairs.length + AtomPairBytes*totalNumAtomPairs);
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
			offset += ResPairBytes + AtomPairBytes*resPairs[i].info.numAtomPairs;
		}
		
		// put the res pairs and atom pairs
		for (int i=0; i<resPairs.length; i++) {
			ResPair resPair = resPairs[i];
			
			// put the res pair
			databuf.putLong(resPair.info.numAtomPairs);
			databuf.putLong(atomOffsetsByResIndex[resPair.resIndex1]);
			databuf.putLong(atomOffsetsByResIndex[resPair.resIndex2]);
			databuf.putDouble(resPair.weight);
			databuf.putDouble(resPair.offset);
			
			// put the atom pairs
			for (int j=0; j<resPair.info.numAtomPairs; j++) {
				databuf.putLong(resPair.info.flags[j]);
				for (int k=0; k<resPair.info.numPrecomputedPerAtomPair; k++) {
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
		func.blockThreads = 1024; // TODO: use func calculator
		func.sharedMemCalc = (int blockThreads) -> blockThreads*Double.BYTES;
	}
	
	@Override
	public void cleanup() {
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
		for (Residue res : ress) {
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
