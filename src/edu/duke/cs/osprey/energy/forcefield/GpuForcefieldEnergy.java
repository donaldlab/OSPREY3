package edu.duke.cs.osprey.energy.forcefield;

import java.io.IOException;
import java.nio.DoubleBuffer;

import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.gpu.kernels.ForceFieldKernel;

public class GpuForcefieldEnergy implements EnergyFunction {
	
	private static final long serialVersionUID = -9142317985561910731L;
	
	private BigForcefieldEnergy ffenergy;
	private ForceFieldKernel.Bound kernel;
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions)
	throws IOException {
		
		ffenergy = new BigForcefieldEnergy(ffparams, interactions, true);
		
		// prep the kernel, upload precomputed data
		kernel = new ForceFieldKernel().bind();
	}
	
	public BigForcefieldEnergy getForcefieldEnergy() {
		return ffenergy;
	}
	
	public int getGpuBytesNeeded() {
		return kernel.getGpuBytesNeeded();
	}
	
	public void initGpu() {
		kernel.setForcefield(ffenergy);
		kernel.uploadStaticAsync();
		kernel.waitForGpu();
	}
	
	@Override
	public double getEnergy() {
		
		// upload coords
		ffenergy.updateCoords();
		kernel.uploadCoordsAsync();
		
		// run the kernel
		kernel.runAsync();
		
		// read the results
		DoubleBuffer out = kernel.downloadEnergiesSync();
		
		// do the last bit of the energy sum on the cpu
		// add one element per work group on the gpu
		// typically, it's a factor of 1024 less than the number of atom pairs
		out.rewind();
		double energy = ffenergy.getInternalSolvationEnergy();
		while (out.hasRemaining()) {
			energy += out.get();
		}
		return energy;
	}
	
	public void cleanup() {
		kernel.cleanup();
	}
}
