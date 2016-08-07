package edu.duke.cs.osprey.energy.forcefield;

import java.io.IOException;
import java.nio.DoubleBuffer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLEvent;
import com.jogamp.opencl.CLEvent.ProfilingCommand;
import com.jogamp.opencl.CLEventList;
import com.jogamp.opencl.CLException;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.gpu.kernels.ForceFieldKernel;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.TimeFormatter;

public class GpuForcefieldEnergy implements EnergyFunction.DecomposableByDof, EnergyFunction.NeedsCleanup {
	
	private static final long serialVersionUID = -9142317985561910731L;
	
	private ForcefieldParams ffparams;
	private ForcefieldInteractions interactions;
	private BigForcefieldEnergy ffenergy;
	private ForceFieldKernel.Bound kernel;
	private Map<List<DegreeOfFreedom>,List<GpuForcefieldEnergy>> decomposedEfuncs;
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions, CLCommandQueue queue)
	throws IOException {
		
		this.ffparams = ffparams;
		this.interactions = interactions;
		
		kernel = new ForceFieldKernel().bind(queue);
		makeForcefield();
		decomposedEfuncs = new HashMap<>();
	}
	
	private void makeForcefield() {
		
		// prep the kernel, upload precomputed data
		ffenergy = new BigForcefieldEnergy(ffparams, interactions, true);
		kernel.setForcefield(ffenergy);
		kernel.uploadStaticAsync();
		kernel.waitForGpu();
	}
	
	public BigForcefieldEnergy getForcefieldEnergy() {
		return ffenergy;
	}
	
	public ForceFieldKernel.Bound getKernel() {
		return kernel;
	}
	
	public void startProfile() {
		kernel.initProfilingEvents();
	}
	
	public String dumpProfile() {
		StringBuilder buf = new StringBuilder();
		CLEventList events = kernel.getProfilingEvents();
		for (CLEvent event : events) {
			try {
				long startNs = event.getProfilingInfo(ProfilingCommand.START);
				long endNs = event.getProfilingInfo(ProfilingCommand.END);
				buf.append(String.format(
					"%s %s\n",
					event.getType(),
					TimeFormatter.format(endNs - startNs, TimeUnit.MICROSECONDS)
				));
			} catch (CLException.CLProfilingInfoNotAvailableException ex) {
				buf.append(String.format("%s (unknown timing)\n", event.getType()));
			}
		}
		kernel.clearProfilingEvents();
		return buf.toString();
	}
	
	@Override
	public double getEnergy() {
		return getEnergySync();
	}
	
	public double getEnergySync() {
		
		// look for residue template changes and rebuild forcefield if needed
		boolean isChanged = false;
		for (AtomGroup[] pair : interactions) {
			for (AtomGroup group : pair) {
				if (group.hasChemicalChange()) {
					isChanged = true;
				}
				group.ackChemicalChange();
			}
		}
		if (isChanged) {
			makeForcefield();
		}
		
		// upload coords
		ffenergy.updateCoords();
		kernel.uploadCoordsAsync();
		
		// run the kernel
		kernel.runAsync();
		
		// read the results
		return sumEnergy(kernel.downloadEnergiesSync());
	}
	
	private double sumEnergy(DoubleBuffer buf) {
		
		// do the last bit of the energy sum on the cpu
		// add one element per work group on the gpu
		// typically, it's a factor of groupSize less than the number of atom pairs
		buf.rewind();
		double energy = ffenergy.getInternalSolvationEnergy();
		while (buf.hasRemaining()) {
			energy += buf.get();
		}
		return energy;
	}
	
	@Override
	public void cleanup() {
		
		kernel.cleanup();
		
		for (List<GpuForcefieldEnergy> efuncs : decomposedEfuncs.values()) {
			for (GpuForcefieldEnergy efunc : efuncs) {
				efunc.cleanup();
			}
		}
	}

	@Override
	public List<EnergyFunction> decomposeByDof(Molecule m, List<DegreeOfFreedom> dofs) {
		
		// check the cache first
		List<GpuForcefieldEnergy> efuncs = decomposedEfuncs.get(dofs);
		if (efuncs == null) {
			
			// cache miss, do the decomposition
			efuncs = new ArrayList<>();
			
			for (DegreeOfFreedom dof : dofs) {

				Residue res = dof.getResidue();
				if (res == null) {
					
					// when there's no residue at the dof, then use the whole efunc
					efuncs.add(this);
					
				} else {
					
					// otherwise, make an efunc for only that residue
					try {
						
						efuncs.add(new GpuForcefieldEnergy(
							ffenergy.getParams(),
							interactions.makeSubsetByResidue(res),
							kernel.getQueue()
						));
						
					} catch (IOException ex) {
						
						// couldn't init gpu kernel for some reason, bail hard
						throw new Error("can't init gpu kernel for decomposed energy function", ex);
					}
				}
			}
			
			// update the cache
			decomposedEfuncs.put(dofs, efuncs);
		}
		
		return new ArrayList<>(efuncs);
	}
}
