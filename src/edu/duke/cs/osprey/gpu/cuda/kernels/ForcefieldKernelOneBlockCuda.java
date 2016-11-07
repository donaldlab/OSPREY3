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
import jcuda.Pointer;

public class ForcefieldKernelOneBlockCuda extends Kernel implements ForcefieldKernel {
	
	private CUBuffer<DoubleBuffer> coords;
	private CUBuffer<IntBuffer> atomFlags;
	private CUBuffer<DoubleBuffer> precomputed;
	private CUBuffer<IntBuffer> subsetTable;
	private CUBuffer<DoubleBuffer> energies;
	
	private CUBuffer<ByteBuffer> args;
	
	private Pointer pKernelArgs;
	
	private int blockThreads;
	private int maxNumBlocks;
	private int numBlocks;
	
	private BigForcefieldEnergy ffenergy;
	private BigForcefieldEnergy.Subset subset;
	
	public ForcefieldKernelOneBlockCuda(GpuStream stream, BigForcefieldEnergy ffenergy)
	throws IOException {
		super(stream, "forcefieldOneBlock", "calcEnergy");
		
		// OPTIMIZATION: 1024 threads seems to work best on modern cards, including GTX 1070 and Tesla K80
		blockThreads = 512;
		numBlocks = 0; // calculated later by setSubset()
		
		// TODO: make arg
		maxNumBlocks = 4;
		
		this.ffenergy = ffenergy;
		
		// wrap the incoming buffers
		coords = stream.makeBuffer(ffenergy.getCoords());
		atomFlags = stream.makeBuffer(ffenergy.getAtomFlags());
		precomputed = stream.makeBuffer(ffenergy.getPrecomputed());
		
		// make the args buffer
		args = stream.makeByteBuffer(32);
		ByteBuffer argsBuf = args.getHostBuffer();
		argsBuf.rewind();
		argsBuf.putDouble(ffenergy.getCoulombFactor());
		argsBuf.putDouble(ffenergy.getScaledCoulombFactor());
		argsBuf.putDouble(ffenergy.getSolvationCutoff2());
		argsBuf.put((byte)(ffenergy.useDistDependentDielectric() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHElectrostatics() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHVdw() ? 1 : 0));
		
		// make the energy buffer for the full subset
		energies = stream.makeDoubleBuffer(calcNumBlocks(ffenergy.getFullSubset()));
		
		// upload static info
		atomFlags.uploadAsync();
		precomputed.uploadAsync();
		args.uploadAsync();
		
		// set the subset
		subsetTable = stream.makeIntBuffer(ffenergy.getFullSubset().getNumAtomPairs());
		subset = null;
		setSubset(ffenergy.getFullSubset());
	}
	
	private int calcNumBlocks(BigForcefieldEnergy.Subset subset) {
		return Math.min(maxNumBlocks, divUp(subset.getNumAtomPairs(), blockThreads));
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
		
		// short circuit: don't change things unless we need to
		if (this.subset == subset) {
			return false;
		}
		
		this.subset = subset;
		this.numBlocks = calcNumBlocks(subset);
		
		// upload subset table
		subsetTable.getHostBuffer().clear();
		subset.getSubsetTable().rewind();
		subsetTable.getHostBuffer().put(subset.getSubsetTable());
		subsetTable.getHostBuffer().flip();
		subsetTable.uploadAsync();
		
		// update kernel args
		pKernelArgs = Pointer.to(
			coords.makeDevicePointer(),
			atomFlags.makeDevicePointer(),
			precomputed.makeDevicePointer(),
			args.makeDevicePointer(),
			subsetTable.makeDevicePointer(),
			Pointer.to(new int[] { subset.getNumAtomPairs() }),
			Pointer.to(new int[] { subset.getNum14AtomPairs() }),
			energies.makeDevicePointer()
		);
		
		return true;
	}
	
	@Override
	public void runAsync() {
		int sharedMemSize = blockThreads*Double.BYTES;
		runAsync(numBlocks, blockThreads, sharedMemSize, pKernelArgs);
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
		}
	}
	
	@Override
	public void uploadCoordsAsync() {
		
		// tell the forcefield to gather updated coords
		ffenergy.updateCoords();
		
		coords.uploadAsync();
	}

	@Override
	public double downloadEnergySync() {
		
		DoubleBuffer buf = energies.downloadSync();
		buf.rewind();
		
		// do the final reduction across blocks on the cpu
		double energy = subset.getInternalSolvationEnergy();
		for (int i=0; i<numBlocks; i++) {
			energy += buf.get();
		}
		return energy;
	}

	// TODO: change interface to suck less
	@Override
	public void setForcefield(BigForcefieldEnergy ffenergy) {
		throw new Error("don't do that");
	}
}
