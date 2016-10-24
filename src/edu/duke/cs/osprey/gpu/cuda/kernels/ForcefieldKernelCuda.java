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

public class ForcefieldKernelCuda extends Kernel implements ForcefieldKernel {
	
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
	
	public ForcefieldKernelCuda(Context context)
	throws IOException {
		super(context, "forcefield.cubin", "calc");
		
		// OPTIMIZATION: 128 is empirically fastest on my GeForce 560 Ti at home
		// but we should experiment with the CS dept Teslas too
		blockThreads = 128;
		
		// init defaults
		numBlocks = 0;
		ffenergy = null;
		subset = null;
	}
	
	public CUBuffer<DoubleBuffer> getCoords() {
		return coords;
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
	public void setForcefield(BigForcefieldEnergy ffenergy) {
		
		if (this.ffenergy != null) {
			throw new IllegalStateException("kernel already has a force field");
		}
		
		this.ffenergy = ffenergy;
		
		// allocate the buffers
		coords = new CUBuffer<>(getContext(), ffenergy.getCoords());
		atomFlags = new CUBuffer<>(getContext(), ffenergy.getAtomFlags());
		precomputed = new CUBuffer<>(getContext(), ffenergy.getPrecomputed());
		subsetTable = new CUBuffer<>(getContext(), BufferTools.makeInt(ffenergy.getFullSubset().getNumAtomPairs(), BufferTools.Type.Direct));
		energies = new CUBuffer<>(getContext(), BufferTools.makeDouble(getEnergySize(ffenergy.getFullSubset(), blockThreads), BufferTools.Type.Direct));
		
		// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
		atomFlags.getHostBuffer().rewind();
		precomputed.getHostBuffer().rewind();
		
		// upload static info
		atomFlags.uploadAsync();
		precomputed.uploadAsync();
		
		// make the args buffer
		args = new CUBuffer<>(getContext(), BufferTools.makeByte(36, BufferTools.Type.Direct));
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
		argsBuf.put((byte)1);
		argsBuf.flip();
		
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
	
	@Override
	public BigForcefieldEnergy.Subset getSubset() {
		return subset;
	}
	
	@Override
	public boolean setSubset(BigForcefieldEnergy.Subset subset) {
		checkInit();
		return setSubsetInternal(subset);
	}
	
	private boolean setSubsetInternal(BigForcefieldEnergy.Subset subset) {
		
		// short circuit: don't change things unless we need to
		if (this.subset == subset) {
			return false;
		}
		
		this.subset = subset;
		
		numBlocks = calcNumBlocks(subset.getNumAtomPairs(), blockThreads);
		
		// update kernel args and upload
		ByteBuffer buf = args.getHostBuffer();
		buf.putInt(0, subset.getNumAtomPairs());
		buf.putInt(4, subset.getNum14AtomPairs());
		buf.put(35, (byte)1); // doEnergy = true
		buf.rewind();
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
		checkInit();
		runAsync(numBlocks, blockThreads, blockThreads*Double.BYTES, pKernelArgs);
	}
	
	private static int getEnergySize(BigForcefieldEnergy.Subset subset, int blockThreads) {
		return calcNumBlocks(subset.getNumAtomPairs(), blockThreads);
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
	
	private void checkInit() {
		if (coords == null) {
			throw new IllegalStateException("call setForcefield() before calling anything else");
		}
	}

	@Override
	public void uploadCoordsAsync() {
		checkInit();
		
		// tell the forcefield to gather updated coords
		ffenergy.updateCoords();
		
		coords.uploadAsync();
	}

	@Override
	public double downloadEnergySync() {
		checkInit();
		
		energies.downloadSync();
		DoubleBuffer buf = energies.getHostBuffer();
		
		// do the last bit of the energy sum on the cpu
		// add one element per work group on the gpu
		// typically, it's a factor of groupSize less than the number of atom pairs
		double energy = subset.getInternalSolvationEnergy();
		buf.rewind();
		int n = getEnergySize(subset, blockThreads);
		for (int i=0; i<n; i++) {
			energy += buf.get();
		}
		return energy;
	}
}
