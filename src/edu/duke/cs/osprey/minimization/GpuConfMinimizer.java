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

package edu.duke.cs.osprey.minimization;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public class GpuConfMinimizer extends ConfMinimizer {
	
	public static enum Type {
		
		Cuda {
			
			@Override
			public boolean isSupported() {
				return !edu.duke.cs.osprey.gpu.cuda.Gpus.get().getGpus().isEmpty();
			}
			
			@Override
			public Context makeContext(int numGpus, int streamsPerGpu) {
				
				// use the Cuda GPU energy function, but do CCD on the CPU
				// (the GPU CCD implementation can't handle non-dihedral dofs yet)
				return new Context() {
					
					private GpuStreamPool pool;
					
					{
						pool = new GpuStreamPool(numGpus, streamsPerGpu);
						minimizers = (mof) -> new SimpleCCDMinimizer(mof);
					}
					
					@Override
					public int getNumStreams() {
						return pool.getNumStreams();
					}
					
					@Override
					public EnergyFunction makeEfunc(ForcefieldParams ffparams, ForcefieldInteractions interactions) {
						return new GpuForcefieldEnergy(ffparams, interactions, pool);
					}
					
					@Override
					public void cleanup() {
						pool.cleanup();
					}
				};
			}
		},
		CudaCCD {
			
			@Override
			public boolean isSupported() {
				return Cuda.isSupported();
			}
			
			@Override
			public Context makeContext(int numGpus, int streamsPerGpu) {
				
				// use a CPU energy function, but send it to the Cuda CCD minimizer
				// (which has a built-in GPU energy function)
				return new Context() {
					
					private GpuStreamPool pool;
					
					{
						pool = new GpuStreamPool(numGpus, streamsPerGpu);
						minimizers = (mof) -> new CudaCCDMinimizer(pool, mof);
					}
					
					@Override
					public int getNumStreams() {
						return pool.getNumStreams();
					}
					
					@Override
					public EnergyFunction makeEfunc(ForcefieldParams ffparams, ForcefieldInteractions interactions) {
						return new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
					}
					
					@Override
					public void cleanup() {
						pool.cleanup();
					}
				};
			}
		},
		OpenCL {
			
			@Override
			public boolean isSupported() {
				return !edu.duke.cs.osprey.gpu.opencl.Gpus.get().getGpus().isEmpty();
			}

			@Override
			public Context makeContext(int numGpus, int streamsPerGpu) {
				
				// use the CPU CCD minimizer, with an OpenCL energy function
				return new Context() {
					
					private GpuQueuePool pool;
					
					{
						pool = new GpuQueuePool(numGpus, streamsPerGpu);
						minimizers = (mof) -> new SimpleCCDMinimizer(mof);
					}
					
					@Override
					public int getNumStreams() {
						return pool.getNumQueues();
					}
					
					@Override
					public EnergyFunction makeEfunc(ForcefieldParams ffparams, ForcefieldInteractions interactions) {
						return new GpuForcefieldEnergy(ffparams, interactions, pool);
					}
					
					@Override
					public void cleanup() {
						pool.cleanup();
					}
				};
			}
		};
		
		private static abstract class Context {
			
			public Factory<? extends EnergyFunction,Molecule> efuncs;
			public Factory<? extends Minimizer,MoleculeModifierAndScorer> minimizers;
			
			public Context() {
				efuncs = null;
				minimizers = null;
			}
			
			public abstract int getNumStreams();
			public abstract EnergyFunction makeEfunc(ForcefieldParams ffparams, ForcefieldInteractions interactions);
			public abstract void cleanup();
		}
		
		public abstract boolean isSupported();
		public abstract Context makeContext(int numGpus, int streamsPerGpu);
		
		public static Type pickBest(ConfSpace confSpace) {
			
			// prefer cuda over opencl, when both are available
			// only because our cuda code is much better than the opencl code right now
			if (Cuda.isSupported()) {
				if (hasAllDihedralDofs(confSpace)) {
					return CudaCCD;
				} else {
					return Cuda;
				}
			}
			
			if (OpenCL.isSupported()) {
				return OpenCL;
			}
			
			return null;
		}
		
		public static Type pickBestOrThrow(ConfSpace confSpace) {
			Type type = Type.pickBest(confSpace);
			if (type == null) {
				throw new Error("GPU computation is not supported on this machine. Use CPU computation instead.");
			}
			return type;
		}
	}
	
	private static boolean hasAllDihedralDofs(ConfSpace confSpace) {
		for (int pos=0; pos<confSpace.numPos; pos++) {
			for (RC rc : confSpace.posFlex.get(pos).RCs) {
				for (DegreeOfFreedom dof : rc.DOFs) {
					if (!(dof instanceof FreeDihedral)) {
						return false;
					}
				}
			}
		}
		return true;
	}
	
	public static class Builder {
		
		public final ForcefieldParams ffparams;
		public final Factory<ForcefieldInteractions,Molecule> interactions;
		public final ConfSpace confSpace;
		
		public Type type;
		public int numGpus;
		public int numStreamsPerGpu;
		
		public Builder(ForcefieldParams ffparams, Factory<ForcefieldInteractions,Molecule> interactions, ConfSpace confSpace) {
			
			this.ffparams = ffparams;
			this.interactions = interactions;
			this.confSpace = confSpace;
			
			type = null;
			numGpus = 1;
			numStreamsPerGpu = 1;
		}
		
		public Builder setGpuInfo(Type type, int numGpus, int numStreamsPerGpu) {
			this.type = type;
			this.numGpus = numGpus;
			this.numStreamsPerGpu = numStreamsPerGpu;
			return this;
		}
		
		public GpuConfMinimizer build() {
			return new GpuConfMinimizer(type, numGpus, numStreamsPerGpu, ffparams, interactions, confSpace);
		}
	}
	
	private Type.Context context;
	
	public GpuConfMinimizer(Type type, int numGpus, int streamsPerGpu, ForcefieldParams ffparams, Factory<ForcefieldInteractions,Molecule> interactions, ConfSpace confSpace) {
		
		if (type == null) {
			type = Type.pickBestOrThrow(confSpace);
		}
		
		// make the gpu context
		context = type.makeContext(numGpus, streamsPerGpu);
		
		// make the minimizer
		Factory<? extends EnergyFunction,Molecule> efuncs = new Factory<EnergyFunction,Molecule>() {
			@Override
			public EnergyFunction make(Molecule mol) {
				return context.makeEfunc(ffparams, interactions.make(mol));
			}
		};
		init(context.getNumStreams(), efuncs, context.minimizers, confSpace);
	}
	
	@Override
	public void cleanup() {
		super.cleanup();
		context.cleanup();
	}
}
