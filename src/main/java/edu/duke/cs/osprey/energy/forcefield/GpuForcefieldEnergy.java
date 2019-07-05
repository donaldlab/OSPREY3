/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.energy.forcefield;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.ForcefieldKernel;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.ForcefieldKernelCuda;
import edu.duke.cs.osprey.gpu.opencl.GpuQueue;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.gpu.opencl.kernels.ForcefieldKernelOpenCL;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class GpuForcefieldEnergy implements EnergyFunction.DecomposableByDof, EnergyFunction.NeedsCleanup, EnergyFunction.ExplicitChemicalChanges {
	
	private static final long serialVersionUID = -9142317985561910731L;
	
	private class KernelBuilder {
		
		private int subsetSequenceNumber;
		private ForcefieldKernel kernel;
		
		public KernelBuilder() {
			subsetSequenceNumber = -1;
			kernel = null;
		}
		
		public ForcefieldKernel get(int expectedSequenceNumber) {
			
			// if this kernel doesn't match the sequence, wipe it
			if (kernel != null && expectedSequenceNumber != subsetSequenceNumber) {
				kernel.cleanup();
				kernel = null;
			}
			
			// make a new kernel if needed
			if (kernel == null) {
				try {
					
					if (openclQueue != null) {
						kernel = new ForcefieldKernelOpenCL(openclQueue, ffenergy);
					} else if (cudaStream != null) {
						kernel = new ForcefieldKernelCuda(cudaStream, ffenergy);
					} else {
						throw new Error("bad gpu queue/context configuration, this is a bug");
					}
					
					// update the sequence
					subsetSequenceNumber = kernel.getForcefield().getFullSubset().handleChemicalChanges();
					
				} catch (IOException ex) {
					
					// if we can't find the gpu kernel source, that's something a programmer needs to fix
					throw new Error("can't initialize gpu kernel", ex);
				}
			}
			
			return kernel;
		}
		
		public void cleanup() {
			if (kernel != null) {
				kernel.cleanup();
				kernel = null;
			}
		}
		
		@Override
		protected void finalize()
		throws Throwable {
			try {
				if (kernel != null) {
					System.err.println("WARNING: " + getClass().getName() + " was garbage collected, but not cleaned up. Attempting cleanup now");
					cleanup();
				}
			} finally {
				super.finalize();
			}
		}
	}
	
	private BigForcefieldEnergy ffenergy;
	private BigForcefieldEnergy.Subset ffsubset;
	private GpuQueuePool openclQueuePool;
	private GpuQueue openclQueue;
	private GpuStreamPool cudaStreamPool;
	private GpuStream cudaStream;
	private KernelBuilder kernelBuilder;
	private Map<Residue,GpuForcefieldEnergy> efuncCache;
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions, GpuQueuePool queuePool) {
		this.ffenergy = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
		this.ffsubset = ffenergy.getFullSubset();
		this.openclQueuePool = queuePool;
		this.openclQueue = queuePool.checkout();
		this.cudaStreamPool = null;
		this.cudaStream = null;
		this.kernelBuilder = new KernelBuilder();
		this.efuncCache = new HashMap<>();
	}
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions, GpuStreamPool streamPool) {
		this.ffenergy = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
		this.ffsubset = ffenergy.getFullSubset();
		this.openclQueuePool = null;
		this.openclQueue = null;
		this.cudaStreamPool = streamPool;
		this.cudaStream = streamPool.checkout();
		this.kernelBuilder = new KernelBuilder();
		this.efuncCache = new HashMap<>();
	}
	
	public GpuForcefieldEnergy(GpuForcefieldEnergy parent, ForcefieldInteractions interactions) {
		this.ffenergy = parent.ffenergy;
		this.ffsubset = ffenergy.new Subset(interactions);
		if (parent.openclQueue != null) {
			this.openclQueuePool = null;
			this.openclQueue = parent.openclQueue;
			this.cudaStreamPool = null;
			this.cudaStream = null;
		} else {
			this.openclQueuePool = null;
			this.openclQueue = null;
			this.cudaStreamPool = null;
			this.cudaStream = parent.cudaStream;
		}
		this.kernelBuilder = parent.kernelBuilder;
		this.efuncCache = null;
	}
	
	public ForcefieldKernel getKernel() {
		return kernelBuilder.get(handleChemicalChanges());
	}
	
	public BigForcefieldEnergy.Subset getSubset() {
		return ffsubset;
	}
	
	@Override
	public int handleChemicalChanges() {
		return ffsubset.handleChemicalChanges();
	}
	
	@Override
	public double getEnergy() {
		
		// check for broken confs
		if (ffsubset.isBroken()) {
			return Double.POSITIVE_INFINITY;
		}
		
		// upload data
		ForcefieldKernel kernel = getKernel();
		kernel.setSubset(getSubset());
		kernel.uploadCoordsAsync();
		
		// compute the energies
		kernel.runAsync();
		
		// read the results
		return kernel.downloadEnergySync();
	}
	
	@Override
	public void clean() {
		
		kernelBuilder.cleanup();
		kernelBuilder = null;
		
		if (openclQueuePool != null) {
			openclQueuePool.release(openclQueue);
		}
		
		if (cudaStreamPool != null) {
			cudaStreamPool.release(cudaStream);
		}
		
		if (efuncCache != null) {
			for (GpuForcefieldEnergy efunc : efuncCache.values()) {
				efunc.clean();
			}
			efuncCache.clear();
		}
	}
	
	@Override
	protected void finalize()
	throws Throwable {
		try {
			if (kernelBuilder != null) {
				System.err.println("WARNING: " + getClass().getName() + " was garbage collected, but not cleaned up. Attempting cleanup now");
				clean();
			}
		} finally {
			super.finalize();
		}
	}

	@Override
	public List<EnergyFunction> decomposeByDof(Molecule m, List<DegreeOfFreedom> dofs) {
		
		List<EnergyFunction> efuncs = new ArrayList<>();
		
		for (DegreeOfFreedom dof : dofs) {

			Residue res = dof.getResidue();
			if (res == null) {
				
				// when there's no residue at the dof, then use the whole efunc
				efuncs.add(this);
				
			} else {
				
				// otherwise, make an efunc for only that residue
				// but share efuncs between dofs in the same residue
				GpuForcefieldEnergy efunc = efuncCache.get(res);
				if (efunc == null) {
					efunc = new GpuForcefieldEnergy(this, ffenergy.getInteractions().makeSubsetByResidue(res));
					efuncCache.put(res, efunc);
				}
				efuncs.add(efunc);
			}
		}
		
		return efuncs;
	}
}
