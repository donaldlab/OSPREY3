package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.ForcefieldKernel;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.Context;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import jcuda.Pointer;

public class ForcefieldKernelOneBlockCuda extends Kernel implements ForcefieldKernel {
	
	private CUBuffer<DoubleBuffer> coords;
	private CUBuffer<IntBuffer> atomFlags;
	private CUBuffer<DoubleBuffer> precomputed;
	private CUBuffer<IntBuffer> subsetTable;
	private CUBuffer<DoubleBuffer> energies;
	
	private CUBuffer<ByteBuffer> args;
	
	private int blockThreads;
	private int numBlocks;
	private Pointer pKernelArgs;
	
	private BigForcefieldEnergy ffenergy;
	private BigForcefieldEnergy.Subset subset;
	
	public ForcefieldKernelOneBlockCuda(Context context, BigForcefieldEnergy ffenergy)
	throws IOException {
		super(context, "forcefieldOneBlock", "calcEnergy");
		
		// OPTIMIZATION: 640 is empirically fastest on my GeForce 560 Ti at home
		// but we should experiment with the CS dept Teslas too
		blockThreads = 512 + 128;
		numBlocks = 1;
		
		this.ffenergy = ffenergy;
		
		// allocate the buffers
		coords = new CUBuffer<>(getContext(), ffenergy.getCoords());
		atomFlags = new CUBuffer<>(getContext(), ffenergy.getAtomFlags());
		precomputed = new CUBuffer<>(getContext(), ffenergy.getPrecomputed());
		subsetTable = new CUBuffer<>(getContext(), BufferTools.makeInt(ffenergy.getFullSubset().getNumAtomPairs(), BufferTools.Type.Direct));
		energies = new CUBuffer<>(getContext(), BufferTools.makeDouble(1, BufferTools.Type.Direct));
		
		// upload static info
		atomFlags.uploadAsync();
		precomputed.uploadAsync();
		
		// make the args buffer
		args = new CUBuffer<>(getContext(), BufferTools.makeByte(40, BufferTools.Type.Direct));
		ByteBuffer argsBuf = args.getHostBuffer();
		argsBuf.rewind();
		argsBuf.putInt(0); // set by setSubsetInternal()
		argsBuf.putInt(0); // 
		argsBuf.putDouble(ffenergy.getCoulombFactor());
		argsBuf.putDouble(ffenergy.getScaledCoulombFactor());
		argsBuf.putDouble(ffenergy.getSolvationCutoff2());
		argsBuf.put((byte)(ffenergy.useDistDependentDielectric() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHElectrostatics() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHVdw() ? 1 : 0));
		
		// set the subset
		// NOTE: setting the subset uploads the args too
		subset = null;
		setSubsetInternal(ffenergy.getFullSubset());
		
		pKernelArgs = Pointer.to(
			coords.makeDevicePointer(),
			atomFlags.makeDevicePointer(),
			precomputed.makeDevicePointer(),
			subsetTable.makeDevicePointer(),
			args.makeDevicePointer(),
			energies.makeDevicePointer()
		);
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
		
		// update kernel args and upload
		ByteBuffer buf = args.getHostBuffer();
		buf.putInt(0, subset.getNumAtomPairs());
		buf.putInt(4, subset.getNum14AtomPairs());
		args.uploadAsync();
		
		// upload subset table
		subsetTable.getHostBuffer().clear();
		subset.getSubsetTable().rewind();
		subsetTable.getHostBuffer().put(subset.getSubsetTable());
		subsetTable.getHostBuffer().flip();
		subsetTable.uploadAsync();
		
		return true;
	}
	
	@Override
	public void runAsync() {
		runAsync(numBlocks, blockThreads, blockThreads*Double.BYTES, pKernelArgs);
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
		return energies.downloadSync().get(0) + subset.getInternalSolvationEnergy();
	}

	// TODO: change interface to suck less
	@Override
	public void setForcefield(BigForcefieldEnergy ffenergy) {
		throw new Error("don't do that");
	}
}
