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
		
		public Bound(Kernel<ForceFieldKernel.Bound> kernel, Gpu gpu) {
			super(kernel, gpu);
		}
		
		public void setForcefield(BigForcefieldEnergy ff) {
			
			int workSize = roundUpWorkSize(ff.getNumAtomPairs());
			setWorkSize(workSize);
			
			CLContext context = getGpu().getDevice().getContext();
			this.coords = context.createBuffer(ff.getCoords(), CLMemory.Mem.READ_WRITE);
			this.atomFlags = context.createBuffer(ff.getAtomFlags(), CLMemory.Mem.READ_ONLY);
			this.precomputed = context.createBuffer(ff.getPrecomputed(), CLMemory.Mem.READ_ONLY);
			this.energies = context.createDoubleBuffer(workSize, CLMemory.Mem.WRITE_ONLY);
			getKernel().getCLKernel()
				.putArg(this.coords)
				.putArg(this.atomFlags)
				.putArg(this.precomputed)
				.putArg(this.energies)
				.putArg(ff.getNumAtomPairs())
				.putArg(ff.getNum14AtomPairs())
				.putArg(ff.getCoulombFactor())
				.putArg(ff.getScaledCoulombFactor())
				.putArg(ff.getSolvationCutoff2())
				.putArg(ff.getSolvationScale())
				// opencl kernels don't support boolean args, so encode as int
				// but bitpack them to save on registers (we're really low on registers in the kernel!)
				.putArg(
					(ff.useDistDependentDielectric() ? 1 : 0)
					| (ff.useHElectrostatics() ? 1 : 0) << 1
					| (ff.useHVdw() ? 1 : 0) << 2
				);
		}
		
		public void uploadStaticAsync() {
			// NOTE: rewind the buffers before uploading, otherwise we get garbage on the gpu
			atomFlags.getBuffer().rewind();
			uploadBufferAsync(atomFlags);
			precomputed.getBuffer().rewind();
			uploadBufferAsync(precomputed);
		}
		
		public void uploadCoordsAsync() {
			coords.getBuffer().rewind();
			uploadBufferAsync(coords);
		}

		public DoubleBuffer downloadEnergiesSync() {
			downloadBufferSync(energies);
			return energies.getBuffer();
		}
	}
}
