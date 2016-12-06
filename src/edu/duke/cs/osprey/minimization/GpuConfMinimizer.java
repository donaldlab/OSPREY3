package edu.duke.cs.osprey.minimization;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public class GpuConfMinimizer extends SpecializedConfMinimizer {
	
	public static enum Type {
		
		Cuda {
			
			@Override
			public boolean isSupported() {
				return !edu.duke.cs.osprey.gpu.cuda.Gpus.get().getGpus().isEmpty();
			}
			
			private boolean hasAllDihedralDofs(ConfSpace confSpace) {
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

			@Override
			public Context makeContext(int numGpus, int streamsPerGpu, ConfSpace confSpace) {
				
				if (hasAllDihedralDofs(confSpace)) {
				
					// use a CPU energy function, but send it to the Cuda CCD minimizer (which has a built-in GPU energy function)
					return new Context() {
						
						private GpuStreamPool pool;
						
						{
							pool = new GpuStreamPool(numGpus, streamsPerGpu);
							minimizers = new Factory<CudaCCDMinimizer,MoleculeModifierAndScorer>() {
								@Override
								public CudaCCDMinimizer make(MoleculeModifierAndScorer mof) {
									CudaCCDMinimizer minimizer = new CudaCCDMinimizer(pool);
									minimizer.init(mof);
									return minimizer;
								}
							};
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
					
				} else {
					
					// use the cuda GPU energy function, but do CCD on the cpu
					// the GPU CCD implementation can't handle non-dihedral dofs yet
					return new Context() {
						
						private GpuStreamPool pool;
						
						{
							pool = new GpuStreamPool(numGpus, streamsPerGpu);
							minimizers = new Factory<SimpleCCDMinimizer,MoleculeModifierAndScorer>() {
								@Override
								public SimpleCCDMinimizer make(MoleculeModifierAndScorer mof) {
									SimpleCCDMinimizer minimizer = new SimpleCCDMinimizer();
									minimizer.init(mof);
									return minimizer;
								}
							};
						}
						
						@Override
						public int getNumStreams() {
							return pool.getNumStreams();
						}
						
						@Override
						public EnergyFunction makeEfunc(ForcefieldParams ffparams, ForcefieldInteractions interactions) {
							return new EnergyFunctionGenerator(ffparams, Double.POSITIVE_INFINITY, false).interactionEnergy(interactions);
						}
						
						@Override
						public void cleanup() {
							pool.cleanup();
						}
					};
				}
			}
		},
		OpenCL {
			
			@Override
			public boolean isSupported() {
				return !edu.duke.cs.osprey.gpu.opencl.Gpus.get().getGpus().isEmpty();
			}

			@Override
			public Context makeContext(int numGpus, int streamsPerGpu, ConfSpace confSpace) {
				
				// use the CPU CCD minimizer, with an OpenCL energy function
				return new Context() {
					
					private GpuQueuePool pool;
					
					{
						pool = new GpuQueuePool(numGpus, streamsPerGpu);
						minimizers = new Factory<SimpleCCDMinimizer,MoleculeModifierAndScorer>() {
							@Override
							public SimpleCCDMinimizer make(MoleculeModifierAndScorer mof) {
								SimpleCCDMinimizer minimizer = new SimpleCCDMinimizer();
								minimizer.init(mof);
								return minimizer;
							}
						};
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
		public abstract Context makeContext(int numGpus, int streamsPerGpu, ConfSpace confSpace);
		
		public static Type pickBest() {
			
			// prefer cuda over opencl, when both are available
			// only because our cuda code is much better than the opencl code right now
			if (Cuda.isSupported()) {
				return Cuda;
			}
			
			if (OpenCL.isSupported()) {
				return OpenCL;
			}
			
			return null;
		}
		
		public static Type pickBestOrThrow() {
			Type type = Type.pickBest();
			if (type == null) {
				throw new Error("GPU computation is not supported on this machine. Use CPU computation instead.");
			}
			return type;
		}
	}
	
	private Type.Context context;
	
	public GpuConfMinimizer(Type type, int numGpus, int streamsPerGpu, ForcefieldParams ffparams, Factory<ForcefieldInteractions,Molecule> interactions, ConfSpace confSpace) {
		
		if (type == null) {
			type = Type.pickBestOrThrow();
		}
		
		// make the gpu context
		context = type.makeContext(numGpus, streamsPerGpu, confSpace);
		
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
