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
import edu.duke.cs.osprey.gpu.Kernel;

public class ForceFieldKernel extends Kernel<ForceFieldKernel.Bound> {
	
	// FYI, useful OpenCL tutorial:
	// http://www.cc.gatech.edu/~vetter/keeneland/tutorial-2011-04-14/06-intro_to_opencl.pdf
	
	// also, reductions on GPUs are tricky. here's a reductions tutorial:
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/
	
	public ForceFieldKernel()
	throws IOException {
		super("forcefield.cl", "calc");
	}
	
	@Override
	public Bound bind(Gpu gpu) {
		return new Bound(this, gpu);
	}
	
	public static class Bound extends BoundKernel<Bound> {
		
		private CLBuffer<DoubleBuffer> coords;
		private CLBuffer<IntBuffer> atomFlags;
		private CLBuffer<DoubleBuffer> precomputed;
		private CLBuffer<DoubleBuffer> energies;
		
		private int workSize;
		private int groupSize;
		
		public Bound(Kernel<ForceFieldKernel.Bound> kernel, Gpu gpu) {
			super(kernel, gpu);
		}
		
		public void setForcefield(BigForcefieldEnergy ffenergy) {
			
			groupSize = getMaxGroupSize();
			workSize = roundUpWorkSize(ffenergy.getNumAtomPairs(), groupSize);
			
			CLContext context = getGpu().getDevice().getContext();
			coords = context.createBuffer(ffenergy.getCoords(), CLMemory.Mem.READ_ONLY);
			atomFlags = context.createBuffer(ffenergy.getAtomFlags(), CLMemory.Mem.READ_ONLY);
			precomputed = context.createBuffer(ffenergy.getPrecomputed(), CLMemory.Mem.READ_ONLY);
			
			energies = context.createDoubleBuffer(workSize/groupSize, CLMemory.Mem.WRITE_ONLY);
			
			getKernel().getCLKernel()
				.setArg(0, coords)
				.setArg(1, atomFlags)
				.setArg(2, precomputed)
				.setArg(3, energies)
				.setArg(4, ffenergy.getNumAtomPairs())
				.setArg(5, ffenergy.getNum14AtomPairs())
				.setArg(6, ffenergy.getCoulombFactor())
				.setArg(7, ffenergy.getScaledCoulombFactor())
				.setArg(8, ffenergy.getSolvationCutoff2())
				// opencl kernels don't support boolean args, so encode as int
				// but bitpack them to save on registers (we're really low on registers in the kernel!)
				.setArg(9,
					(ffenergy.useDistDependentDielectric() ? 1 : 0)
					| (ffenergy.useHElectrostatics() ? 1 : 0) << 1
					| (ffenergy.useHVdw() ? 1 : 0) << 2
				)
				// allocate the gpu local memory
				.setNullArg(10, groupSize*Double.BYTES);
		}
		
		public void uploadStaticAsync() {
			
			// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
			atomFlags.getBuffer().rewind();
			precomputed.getBuffer().rewind();
			
			uploadBufferAsync(atomFlags);
			uploadBufferAsync(precomputed);
		}
		
		public void uploadCoordsAsync() {
			
			// IMPORTANT: rewind the buffers before uploading, otherwise we get garbage on the gpu
			coords.getBuffer().rewind();
			
			uploadBufferAsync(coords);
		}

		public void runAsync() {
			runAsync(workSize, groupSize);
		}
		
		public DoubleBuffer downloadEnergiesSync() {
			
			// IMPORTANT!! rewind the output buffer before downloading energies
			// otherwise we get weird segfaults in nvidia's opencl driver that are next to impossible to diagnose!
			energies.getBuffer().rewind();
			
			downloadBufferSync(energies);
			return energies.getBuffer();
		}
		
		public int getGpuBytesNeeded() {
			return this.coords.getCLCapacity()
				+ this.atomFlags.getCLCapacity()
				+ this.precomputed.getCLCapacity()
				+ workSize/groupSize*Double.BYTES
				+ groupSize*Double.BYTES;
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
	}
}
