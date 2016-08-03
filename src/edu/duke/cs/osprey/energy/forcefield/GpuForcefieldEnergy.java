package edu.duke.cs.osprey.energy.forcefield;

import java.io.IOException;
import java.nio.DoubleBuffer;
import java.util.concurrent.TimeUnit;

import com.jogamp.opencl.CLEvent;
import com.jogamp.opencl.CLEvent.ProfilingCommand;
import com.jogamp.opencl.CLEventList;

import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.gpu.kernels.ForceFieldKernel;
import edu.duke.cs.osprey.tools.TimeFormatter;

public class GpuForcefieldEnergy implements EnergyFunction {
	
	private static final long serialVersionUID = -9142317985561910731L;
	
	private BigForcefieldEnergy ffenergy;
	private ForceFieldKernel.Bound kernel;
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions)
	throws IOException {
		this(ffparams, interactions, false);
	}
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions, boolean useProfiling)
	throws IOException {
		
		ffenergy = new BigForcefieldEnergy(ffparams, interactions, true);
		
		// prep the kernel, upload precomputed data
		kernel = new ForceFieldKernel().bind(useProfiling);
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
	
	public void startProfile() {
		kernel.initProfilingEvents();
	}
	
	public String dumpProfile() {
		StringBuilder buf = new StringBuilder();
		CLEventList events = kernel.getProfilingEvents();
		for (CLEvent event : events) {
			long startNs = event.getProfilingInfo(ProfilingCommand.START);
			long endNs = event.getProfilingInfo(ProfilingCommand.END);
			buf.append(String.format(
				"%s %s\n",
				event.getType(),
				TimeFormatter.format(endNs - startNs, TimeUnit.MICROSECONDS)
			));
		}
		kernel.clearProfilingEvents();
		return buf.toString();
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
