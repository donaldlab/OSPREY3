package edu.duke.cs.osprey.gpu.opencl.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLMemory;

import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.gpu.ForcefieldKernel;
import edu.duke.cs.osprey.gpu.opencl.GpuQueue;
import edu.duke.cs.osprey.gpu.opencl.Kernel;

public class ForcefieldKernelOpenCL extends Kernel implements ForcefieldKernel {
	
	// FYI, useful OpenCL tutorial:
	// http://www.cc.gatech.edu/~vetter/keeneland/tutorial-2011-04-14/06-intro_to_opencl.pdf
	
	// also, reductions on GPUs are tricky. here's a reductions tutorial:
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/
	
	private CLBuffer<DoubleBuffer> coords;
	private CLBuffer<IntBuffer> atomFlags;
	private CLBuffer<DoubleBuffer> precomputed;
	private CLBuffer<IntBuffer> subsetTable;
	private CLBuffer<DoubleBuffer> energies;
	
	private CLBuffer<ByteBuffer> args;
	
	private int workSize;
	private int groupSize;
	
	private BigForcefieldEnergy ffenergy;
	private BigForcefieldEnergy.Subset subset;

	public ForcefieldKernelOpenCL(GpuQueue queue, BigForcefieldEnergy ffenergy)
	throws IOException {
		super(queue, "forcefield.cl", "calc");
		
		/* OPTIMIZATION: this kernel uses lots and lots of registers, so maxing out the work group size is sub-optimal
			using a smaller group size works noticeably better!
			empirically, using 1/4 the max seems to make the difference between 11x speedups and 15x speedups
		
			times in us for (131k atom pairs, 11k atom pairs, 200 atom pairs) on 1000th run, 6000th run, and 100000th run respectively
				1024: 369, 44, 18 (max on my hardware)
				256: 272, 33, 15
				128: 257, 31, 14
				64: 293, 33, 13
				32: 411, 44, 12
			
			looks to have slight dependency on num atom pairs, but 128 looks like a good compromise for all sizes
		*/
		groupSize = 128;
		
		// init defaults
		workSize = 0;
		
		this.ffenergy = ffenergy;
		
		// allocate the buffers
		CLContext context = getQueue().getCLQueue().getContext();
		coords = context.createBuffer(ffenergy.getCoords(), CLMemory.Mem.READ_WRITE);
		atomFlags = context.createBuffer(ffenergy.getAtomFlags(), CLMemory.Mem.READ_ONLY);
		precomputed = context.createBuffer(ffenergy.getPrecomputed(), CLMemory.Mem.READ_ONLY);
		subsetTable = context.createIntBuffer(ffenergy.getFullSubset().getNumAtomPairs(), CLMemory.Mem.READ_ONLY);
		energies = context.createDoubleBuffer(getEnergySize(ffenergy.getFullSubset()), CLMemory.Mem.WRITE_ONLY);
		
		// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
		atomFlags.getBuffer().rewind();
		precomputed.getBuffer().rewind();
		
		// upload static info
		uploadBufferAsync(atomFlags);
		uploadBufferAsync(precomputed);
		
		// make the args buffer
		args = context.createByteBuffer(40, CLMemory.Mem.READ_ONLY);
		ByteBuffer argsBuf = args.getBuffer();
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
		
		getCLKernel()
			.setArg(0, coords)
			.setArg(1, atomFlags)
			.setArg(2, precomputed)
			.setArg(3, subsetTable)
			.setArg(4, energies)
			.setArg(5, args)
			.setNullArg(6, groupSize*Double.BYTES);
	}
	
	public CLBuffer<DoubleBuffer> getCoords() {
		return coords;
	}
	
	public CLBuffer<DoubleBuffer> getEnergies() {
		return energies;
	}
	
	public CLBuffer<ByteBuffer> getArgs() {
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
		
		// make sure the total work size is a multiple of the group size
		workSize = roundUpWorkSize(subset.getNumAtomPairs(), groupSize);
		
		// update kernel args and upload
		ByteBuffer buf = args.getBuffer();
		buf.putInt(0, subset.getNumAtomPairs());
		buf.putInt(4, subset.getNum14AtomPairs());
		buf.put(35, (byte)(useSubset ? 1 : 0));
		buf.rewind();
		uploadBufferAsync(args);
		
		if (useSubset) {
			
			// upload subset table
			subsetTable.getBuffer().clear();
			subset.getSubsetTable().rewind();
			subsetTable.getBuffer().put(subset.getSubsetTable());
			subsetTable.getBuffer().flip();
			uploadPartialBufferAsync(subsetTable);
		}
		
		return true;
	}
	
	public boolean isDoEnergy() {
		downloadBufferSync(args);
		return args.getBuffer().get(35) == 1;
	}
	
	public void setDoEnergy(boolean val) {
		ByteBuffer buf = args.getBuffer();
		buf.put(35, (byte)(val ? 1 : 0));
		buf.rewind();
		uploadBufferAsync(args);
	}
	
	@Override
	public void uploadCoordsAsync() {
		
		// tell the forcefield to gather updated coords
		ffenergy.updateCoords();
	
		// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
		coords.getBuffer().rewind();
		
		uploadBufferAsync(coords);
	}
	
	public DoubleBuffer downloadCoordsSync() {
		coords.getBuffer().rewind();
		downloadBufferSync(coords);
		return coords.getBuffer();
	}

	@Override
	public void runAsync() {
		runAsync(workSize, groupSize);
	}
	
	@Override
	public double downloadEnergySync() {
		
		// IMPORTANT!! rewind the output buffer before downloading energies
		// otherwise we get weird segfaults in nvidia's opencl driver that are next to impossible to diagnose!
		energies.getBuffer().rewind();
		downloadBufferSync(energies);
		
		DoubleBuffer buf = energies.getBuffer();
		
		// do the last bit of the energy sum on the cpu
		// add one element per work group on the gpu
		// typically, it's a factor of groupSize less than the number of atom pairs
		double energy = subset.getInternalSolvationEnergy();
		buf.rewind();
		int n = getEnergySize();
		for (int i=0; i<n; i++) {
			energy += buf.get();
		}
		return energy;
	}
	
	public int getEnergySize() {
		return getEnergySize(subset);
	}
	
	public int getEnergySize(BigForcefieldEnergy.Subset subset) {
		return roundUpWorkSize(subset.getNumAtomPairs(), groupSize)/groupSize;
	}
	
	public int getGpuBytesNeeded() {
		return coords.getCLCapacity()*Double.BYTES
			+ atomFlags.getCLCapacity()*Integer.BYTES
			+ precomputed.getCLCapacity()*Double.BYTES
			+ subsetTable.getCLCapacity()*Integer.BYTES
			+ args.getCLCapacity()
			+ getEnergySize()*Double.BYTES // global kernel memory to save the reduction results
			+ groupSize*Double.BYTES; // local kernel memory for the reduction
	}
	
	@Override
	public void cleanup() {
		
		if (coords != null) {
			
			coords.release();
			atomFlags.release();
			precomputed.release();
			subsetTable.release();
			args.release();
			energies.release();
			
			coords = null;
		}
		
		super.cleanup();
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
