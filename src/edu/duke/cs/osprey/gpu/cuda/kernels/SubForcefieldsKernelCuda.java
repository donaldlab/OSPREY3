package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.structure.Residue;
import jcuda.Pointer;

public class SubForcefieldsKernelCuda extends Kernel {
	
	private CUBuffer<DoubleBuffer> dihedrals;
	private CUBuffer<IntBuffer> dihedralIndices;
	private CUBuffer<IntBuffer> rotatedIndices;
	private CUBuffer<IntBuffer> subsets;
	private CUBuffer<IntBuffer> subsetOffsets;
	private CUBuffer<DoubleBuffer> energies;
	
	private int blockThreads;
	private int numBlocks;
	private Pointer pKernelArgs;
	
	private ForcefieldKernelCuda ffKernel;
	
	public SubForcefieldsKernelCuda(ForcefieldKernelCuda ffKernel, List<FreeDihedral> dofs)
	throws IOException {
		super(ffKernel.getContext(), "subForcefields", "calcEnergies");
		
		this.ffKernel = ffKernel;
		
		// TODO: optimize block threads
		blockThreads = 512;
		numBlocks = 3*dofs.size();
		
		// get the max num rotated indices
		int maxRotatedAtoms = 0;
		for (FreeDihedral dof : dofs) {
			int rotatedAtoms = dof.getResidue().template.getDihedralRotatedAtoms(dof.getDihedralNumber()).size();
			maxRotatedAtoms = Math.max(maxRotatedAtoms, rotatedAtoms);
		}
		
		// get the subset for each unique residue
		Map<Residue,Integer> subsetIndexByResidue = new HashMap<>();
		List<BigForcefieldEnergy.Subset> subsetsList = new ArrayList<>();
		List<Integer> subsetIndexByDof = new ArrayList<>();
		int maxSubsetSize = 0;
		for (FreeDihedral dof : dofs) {
			
			Residue res = dof.getResidue();
			Integer subsetIndex = subsetIndexByResidue.get(res);
			if (subsetIndex == null) {
				
				// make a new subset
				subsetIndex = subsetsList.size();
				ForcefieldInteractions subInteractions = ffKernel.getForcefield().getInteractions().makeSubsetByResidue(res);
				BigForcefieldEnergy.Subset subset = ffKernel.getForcefield().new Subset(subInteractions);
				subsetsList.add(subset);
				subsetIndexByResidue.put(res, subsetIndex);
				
				maxSubsetSize = Math.max(maxSubsetSize, subset.getNumAtomPairs());
			}
			
			subsetIndexByDof.add(subsetIndex);
		}
		
		// allocate the buffers
		dihedrals = new CUBuffer<>(getContext(), BufferTools.makeDouble(numBlocks, BufferTools.Type.Direct));
		dihedralIndices = new CUBuffer<>(getContext(), BufferTools.makeInt(dofs.size()*4, BufferTools.Type.Direct));
		rotatedIndices = new CUBuffer<>(getContext(), BufferTools.makeInt(1 + dofs.size()*(1 + maxRotatedAtoms), BufferTools.Type.Direct));
		subsets = new CUBuffer<>(getContext(), BufferTools.makeInt(subsetsList.size()*(2 + maxSubsetSize), BufferTools.Type.Direct));
		subsetOffsets = new CUBuffer<>(getContext(), BufferTools.makeInt(numBlocks, BufferTools.Type.Direct));
		energies = new CUBuffer<>(getContext(), BufferTools.makeDouble(numBlocks, BufferTools.Type.Direct));
		
		// make the dihedral and rotated indices
		IntBuffer dihedralIndicesBuf = dihedralIndices.getHostBuffer();
		dihedralIndicesBuf.rewind();
		IntBuffer rotatedIndicesBuf = rotatedIndices.getHostBuffer();
		rotatedIndicesBuf.rewind();
		rotatedIndicesBuf.put(maxRotatedAtoms);
		for (FreeDihedral dof : dofs) {
			
			Residue res = dof.getResidue();
			int coordsOffset = ffKernel.getForcefield().getAtomOffset(res);
		
			int[] dihedralIndicesSrc = res.template.getDihedralDefiningAtoms(dof.getDihedralNumber());
			assert (dihedralIndicesSrc.length == 4);
			for (int i=0; i<4; i++) {
				dihedralIndicesBuf.put(dihedralIndicesSrc[i] + coordsOffset);
			}
			
			List<Integer> rotatedIndicesSrc = res.template.getDihedralRotatedAtoms(dof.getDihedralNumber());
			rotatedIndicesBuf.put(rotatedIndicesSrc.size());
			int i;
			for (i=0; i<rotatedIndicesSrc.size(); i++) {
				rotatedIndicesBuf.put(rotatedIndicesSrc.get(i) + coordsOffset);
			}
			for (; i<maxRotatedAtoms; i++) {
				rotatedIndicesBuf.put(-1);
			}
		}
		assert (rotatedIndicesBuf.limit() == rotatedIndicesBuf.capacity());
		
		// make the subsets
		IntBuffer subsetsBuf = subsets.getHostBuffer();
		subsetsBuf.rewind();
		for (BigForcefieldEnergy.Subset subset : subsetsList) {
			subsetsBuf.put(subset.getNumAtomPairs());
			subsetsBuf.put(subset.getNum14AtomPairs());
			subset.getSubsetTable().rewind();
			subsetsBuf.put(subset.getSubsetTable());
			
			int spaceLeft = maxSubsetSize - subset.getNumAtomPairs();
			for (int i=0; i<spaceLeft; i++) {
				subsetsBuf.put(-1);
			}
		}
		assert (subsetsBuf.limit() == subsetsBuf.capacity());
		
		// make the subset offsets
		IntBuffer subsetMapBuf = subsetOffsets.getHostBuffer();
		subsetMapBuf.rewind();
		for (int i=0; i<numBlocks; i++) {
			int subsetIndex = subsetIndexByDof.get(i/3);
			subsetMapBuf.put(subsetIndex*(2 + maxSubsetSize));
		}
		
		// upload static info
		dihedralIndices.uploadAsync();
		rotatedIndices.uploadAsync();
		subsets.uploadAsync();
		subsetOffsets.uploadAsync();
		
		pKernelArgs = Pointer.to(
			ffKernel.getCoords().makeDevicePointer(),
			Pointer.to(new int[] { ffKernel.getCoords().getHostBuffer().capacity() }),
			ffKernel.getAtomFlags().makeDevicePointer(),
			ffKernel.getPrecomputed().makeDevicePointer(),
			ffKernel.getArgs().makeDevicePointer(),
			dihedrals.makeDevicePointer(),
			dihedralIndices.makeDevicePointer(),
			rotatedIndices.makeDevicePointer(),
			subsets.makeDevicePointer(),
			subsetOffsets.makeDevicePointer(),
			energies.makeDevicePointer()
		);
	}
	
	public void calcEnergies(DoubleMatrix1D x, double delta) {
		
		// upload coords
		ffKernel.uploadCoordsAsync();
		
		// upload dihedrals
		DoubleBuffer dihedralsBuf = dihedrals.getHostBuffer();
		dihedralsBuf.rewind();
		for (int d=0; d<x.size(); d++) {
			double xd = x.get(d);
			dihedralsBuf.put(Math.toRadians(xd));
			dihedralsBuf.put(Math.toRadians(xd - delta));
			dihedralsBuf.put(Math.toRadians(xd + delta));
		}
		dihedrals.uploadAsync();
		
		// calc shared memory size
		int sharedMemBytes = 0;
		sharedMemBytes += blockThreads*Double.BYTES; // energy sum
		sharedMemBytes += ffKernel.getCoords().getNumBytes(); // coords copy
		
		runAsync(numBlocks, blockThreads, sharedMemBytes, pKernelArgs);
		
		energies.downloadSync();
		
		/* TEMP: show energies
		DoubleBuffer buf = energies.getHostBuffer();
		buf.rewind();
		for (int d=0; d<x.size(); d++) {
			System.out.println(String.format("\td: %2d, energies: %12.6f   %12.6f   %12.6f", d, buf.get(), buf.get(), buf.get()));
		}
		*/
	}

	public double getFxdmOffset(int d) {
		DoubleBuffer buf = energies.getHostBuffer();
		buf.position(d*3);
		double fxd = buf.get();
		double fxdm = buf.get();
		return fxdm - fxd;
	}

	public double getFxdpOffset(int d) {
		DoubleBuffer buf = energies.getHostBuffer();
		buf.position(d*3);
		double fxd = buf.get();
		buf.get();
		double fxdp = buf.get();
		return fxdp - fxd;
	}
	
	public void cleanup() {
		dihedrals.cleanup();
		dihedralIndices.cleanup();
		rotatedIndices.cleanup();
		subsets.cleanup();
		subsetOffsets.cleanup();
		energies.cleanup();
	}
}
