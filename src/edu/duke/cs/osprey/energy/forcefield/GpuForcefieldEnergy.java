package edu.duke.cs.osprey.energy.forcefield;

import java.io.IOException;
import java.nio.DoubleBuffer;

import edu.duke.cs.osprey.gpu.kernels.ForceFieldKernel;

public class GpuForcefieldEnergy {
	
	private BigForcefieldEnergy ffenergy;
	private ForceFieldKernel.Bound kernel;
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions)
	throws IOException {
		
		ffenergy = new BigForcefieldEnergy(ffparams, interactions, true);
		
		// prep the kernel, upload precomputed data
		kernel = new ForceFieldKernel().bind();
	}
	
	public void initGpu() {
		kernel.setForcefield(ffenergy);
		kernel.uploadStaticAsync();
		kernel.waitForGpu();
	}
	
	public double calculateTotalEnergy() {
		
		// upload coords
		ffenergy.updateCoords();
		kernel.uploadCoordsAsync();
		
		// run the kernel
		kernel.runAsync();
		
		// read the results
		double energy = ffenergy.getInternalSolvationEnergy();
		DoubleBuffer out = kernel.downloadEnergiesSync();
		
		// do the last bit of the energy sum on the cpu
		// add one element per work group on the gpu
		// typically, a factor of 1024 less than the number of atom pairs
		out.rewind();
		while (out.hasRemaining()) {
			energy += out.get();
		}
		return energy;
	}
	
	public void cleanup() {
		kernel.cleanup();
	}
}
