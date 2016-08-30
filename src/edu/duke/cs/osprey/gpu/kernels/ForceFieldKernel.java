package edu.duke.cs.osprey.gpu.kernels;

import java.io.IOException;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLMemory;

import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.gpu.BoundKernel;
import edu.duke.cs.osprey.gpu.Gpu;
import edu.duke.cs.osprey.gpu.GpuQueue;
import edu.duke.cs.osprey.gpu.Kernel;

public class ForceFieldKernel extends Kernel<ForceFieldKernel.Bound> {
	
	// FYI, useful OpenCL tutorial:
	// http://www.cc.gatech.edu/~vetter/keeneland/tutorial-2011-04-14/06-intro_to_opencl.pdf
	
	// also, reductions on GPUs are tricky. here's a reductions tutorial:
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/
	
	public ForceFieldKernel(Gpu gpu)
	throws IOException {
		super(gpu, "forcefield.cl", "calc");
	}
	
	@Override
	public Bound bind(GpuQueue queue) {
		return new Bound(this, queue);
	}
	
	public static class Bound extends BoundKernel<Bound> {
		
		private CLBuffer<DoubleBuffer> coords;
		private CLBuffer<IntBuffer> atomFlags;
		private CLBuffer<DoubleBuffer> precomputed;
		private CLBuffer<IntBuffer> subsetTable;
		private CLBuffer<DoubleBuffer> energies;
		
		private int workSize;
		private int groupSize;
		
		private BigForcefieldEnergy ffenergy;
		private BigForcefieldEnergy.Subset subset;
		
		// NEXTTIME: don't make child kernels, just upload different subset tables when switching energy functions!
		
		public Bound(Kernel<ForceFieldKernel.Bound> kernel, GpuQueue queue) {
			super(kernel, queue);
			
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
			
			// set the subset
			subset = null;
			setSubsetInternal(ffenergy.getFullSubset());
			
			// allocate the energy buffer after setting the full subset
			// this makes sure we have the biggest energy buffer we'll ever need
			energies = context.createDoubleBuffer(getEnergySize(), CLMemory.Mem.WRITE_ONLY);
			
			// set kernel args
			getKernel().getCLKernel()
				.setArg(0, coords)
				.setArg(1, atomFlags)
				.setArg(2, precomputed)
				.setArg(3, subsetTable)
				.setArg(4, energies)
				// args 4 and 5 set by setSubset()
				.setArg(7, ffenergy.getCoulombFactor())
				.setArg(8, ffenergy.getScaledCoulombFactor())
				.setArg(9, ffenergy.getSolvationCutoff2())
				// opencl kernels don't support boolean args, so encode as int
				// but bitpack them to save on registers (we're really low on registers in the kernel!)
				.setArg(10,
					(ffenergy.useDistDependentDielectric() ? 1 : 0)
					| (ffenergy.useHElectrostatics() ? 1 : 0) << 1
					| (ffenergy.useHVdw() ? 1 : 0) << 2
				)
				// allocate the gpu local memory
				.setNullArg(11, groupSize*Double.BYTES);
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
			
			// update kernel args
			getKernel().getCLKernel()
				.setArg(5, subset.getNumAtomPairs())
				.setArg(6, subset.getNum14AtomPairs());
			
			// copy the subset table to the upload buffer
			subsetTable.getBuffer().clear();
			subset.getSubsetTable().rewind();
			subsetTable.getBuffer().put(subset.getSubsetTable());
			
			// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
			subsetTable.getBuffer().rewind();
			
			uploadBufferAsync(subsetTable);
		}
		
		public void uploadStaticAsync() {
			checkInit();
			
			// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
			atomFlags.getBuffer().rewind();
			precomputed.getBuffer().rewind();
			
			uploadBufferAsync(atomFlags);
			uploadBufferAsync(precomputed);
		}
		
		public void uploadCoordsAsync() {
			checkInit();
			
			// tell the forcefield to gather updated coords
			ffenergy.updateCoords();
			
			// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
			coords.getBuffer().rewind();
			
			uploadBufferAsync(coords);
		}

		public void runAsync() {
			checkInit();
			runAsync(workSize, groupSize);
		}
		
		public DoubleBuffer downloadEnergiesSync() {
			checkInit();
			
			// IMPORTANT!! rewind the output buffer before downloading energies
			// otherwise we get weird segfaults in nvidia's opencl driver that are next to impossible to diagnose!
			energies.getBuffer().rewind();
			
			downloadBufferSync(energies);
			return energies.getBuffer();
		}
		
		public int getEnergySize() {
			return workSize/groupSize;
		}
		
		public int getGpuBytesNeeded() {
			return coords.getCLCapacity()*Double.BYTES
				+ atomFlags.getCLCapacity()*Integer.BYTES
				+ precomputed.getCLCapacity()*Double.BYTES
				+ subsetTable.getCLCapacity()*Integer.BYTES
				+ getEnergySize()*Double.BYTES // global kernel memory to save the reduction results
				+ groupSize*Double.BYTES; // local kernel memory for the reduction
		}
		
		@Override
		public void cleanup() {
			coords.release();
			atomFlags.release();
			precomputed.release();
			if (energies != null) {
				energies.release();
			}
		}
		
		private void checkInit() {
			if (coords == null) {
				throw new IllegalStateException("call setForcefield() before calling anything else");
			}
		}
	}
}
