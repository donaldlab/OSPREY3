package edu.duke.cs.osprey.gpu.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLMemory;

import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.gpu.Kernel;
import edu.duke.cs.osprey.gpu.GpuQueue;

public class ForceFieldKernel extends Kernel {
	
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

	public ForceFieldKernel(GpuQueue queue)
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
		ffenergy = null;
		subset = null;
	}
	
	public CLBuffer<DoubleBuffer> getCoords() {
		return coords;
	}
	
	public BigForcefieldEnergy getForcefield() {
		return ffenergy;
	}
	public void setForcefield(BigForcefieldEnergy ffenergy) {
		
		this.ffenergy = ffenergy;
		
		// allocate some of the buffers we need
		CLContext context = getQueue().getCLQueue().getContext();
		coords = context.createBuffer(ffenergy.getCoords(), CLMemory.Mem.READ_ONLY);
		atomFlags = context.createBuffer(ffenergy.getAtomFlags(), CLMemory.Mem.READ_ONLY);
		precomputed = context.createBuffer(ffenergy.getPrecomputed(), CLMemory.Mem.READ_ONLY);
		subsetTable = context.createIntBuffer(ffenergy.getFullSubset().getNumAtomPairs(), CLMemory.Mem.READ_ONLY);
		
		// make the args buffer
		args = context.createByteBuffer(35, CLMemory.Mem.READ_ONLY);
		ByteBuffer argsBuf = args.getBuffer();
		argsBuf.rewind();
		argsBuf.putInt(0); // set by setSubsetInternal()
		argsBuf.putInt(0); // 
		argsBuf.putDouble(ffenergy.getCoulombFactor());
		argsBuf.putDouble(ffenergy.getScaledCoulombFactor());
		argsBuf.putDouble(ffenergy.getSolvationCutoff2());
		argsBuf.put((byte)(ffenergy.useDistDependentDielectric() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHElectrostatics() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHVdw() ? 1 : 0));
		argsBuf.flip();
		
		// set the subset
		subset = null;
		setSubsetInternal(ffenergy.getFullSubset());
		
		// allocate the energy buffer after setting the full subset
		// this makes sure we have the biggest energy buffer we'll ever need
		energies = context.createDoubleBuffer(getEnergySize(), CLMemory.Mem.WRITE_ONLY);
		
		getCLKernel()
			.setArg(0, coords)
			.setArg(1, atomFlags)
			.setArg(2, precomputed)
			.setArg(3, subsetTable)
			.setArg(4, energies)
			.setArg(5, args)
			.setNullArg(6, groupSize*Double.BYTES);
	}
	
	public BigForcefieldEnergy.Subset getSubset() {
		return subset;
	}
	
	public void setSubset(BigForcefieldEnergy.Subset subset) {
		checkInit();
		setSubsetInternal(subset);
	}
	
	private void setSubsetInternal(BigForcefieldEnergy.Subset subset) {
		
		// short circuit: don't change things unless we need to
		if (this.subset == subset) {
			return;
		}
		
		this.subset = subset;
		
		// make sure the total work size is a multiple of the group size
		workSize = roundUpWorkSize(subset.getNumAtomPairs(), groupSize);
		
		// update kernel args and upload
		ByteBuffer argsBuf = args.getBuffer();
		argsBuf.rewind();
		argsBuf.putInt(subset.getNumAtomPairs());
		argsBuf.putInt(subset.getNum14AtomPairs());
		argsBuf.rewind();
		uploadBufferAsync(args);
		
		// upload subset table
		subsetTable.getBuffer().clear();
		subset.getSubsetTable().rewind();
		subsetTable.getBuffer().put(subset.getSubsetTable());
		subsetTable.getBuffer().rewind();
		uploadBufferAsync(subsetTable);
	}
	
	public void uploadStaticAsync() {
		checkInit();
		
		// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
		atomFlags.getBuffer().rewind();
		precomputed.getBuffer().rewind();
		args.getBuffer().rewind();
		
		uploadBufferAsync(atomFlags);
		uploadBufferAsync(precomputed);
		uploadBufferAsync(args);
	}
	
	public void updateAndUploadCoordsAsync() {
		checkInit();
		
		// tell the forcefield to gather updated coords
		ffenergy.updateCoords();
	
		uploadCoordsAsync();
	}
	
	public void uploadCoordsAsync() {
		
		// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
		coords.getBuffer().rewind();
		
		uploadBufferAsync(coords);
	}

	public void runAsync() {
		checkInit();
		runAsync(workSize, groupSize);
	}
	
	public double downloadEnergySync() {
		checkInit();
		
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
		return workSize/groupSize;
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
		coords.release();
		atomFlags.release();
		precomputed.release();
		subsetTable.release();
		args.release();
		energies.release();
	}
	
	private void checkInit() {
		if (coords == null) {
			throw new IllegalStateException("call setForcefield() before calling anything else");
		}
	}
}
