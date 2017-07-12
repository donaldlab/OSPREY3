package edu.duke.cs.osprey.energy;

import java.util.function.Consumer;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.DofTypes;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.ResidueCudaCCDMinimizer;
import edu.duke.cs.osprey.gpu.cuda.kernels.ResidueForcefieldEnergyCuda;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.CudaCCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residues;
import edu.duke.cs.osprey.tools.AutoCleanable;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.UseableBuilder;

/**
 * Computes the energy of a molecule fragment from a conformation space using the desired forcefield parameters
 * and residue interactions.
 * 
 * Residue interactions are specified via {@link ResidueInteractions} instances.
 * Forcefield implementations are chosen via the type argument. If no type is specified, the best forcefield
 * available implementation is automatically chosen based on the parallelism argument.
 * 
 * If a fragment has continuous degrees of freedom, minimization will be performed before forcefield evaluation.
 */
public class EnergyCalculator implements AutoCleanable {
	
	public static class Builder implements UseableBuilder<EnergyCalculator> {
		
		private ForcefieldParams ffparams;
		
		/**
		 * The parallel hardware available for Osprey calculations.
		 * CPU thread pools are available, as well as high-performance GPU kernels.
		 * Just describe your available hardware, and Osprey will try to pick the best implementation.
		 * If you need more control, try the setType() option directly.
		 */
		private Parallelism parallelism = Parallelism.makeCpu(1);
		
		/**
		 * Directly choose the energy calculator implementation.
		 * 
		 * @warning Not all energy calculator implementations may be compatible with your design's
		 * conformation space. Choosing an incompatible implementation can cause Osprey to crash.
		 * 
		 * If you're unsure which implementation to use, just use setParallelism() instead of setting
		 * the type directly.
		 */
		private Type type = null;
		
		
		private DofTypes dofTypes = DofTypes.Any;
		private AtomConnectivity.Builder atomConnectivityBuilder = new AtomConnectivity.Builder();
		private ResPairCache resPairCache;
		
		public Builder(SimpleConfSpace confSpace, ForcefieldParams ffparams) {
			this.ffparams = ffparams;
			this.atomConnectivityBuilder.addTemplates(confSpace);
			this.dofTypes = confSpace.getDofTypes();
		}
		
		public Builder(Residues residues, ForcefieldParams ffparams) {
			this.ffparams = ffparams;
			this.atomConnectivityBuilder.addTemplates(residues);
		}
		
		public Builder setParallelism(Parallelism val) {
			parallelism = val;
			return this;
		}
		
		public Builder setType(Type val) {
			type = val;
			return this;
		}
		
		public Builder setDofTypes(DofTypes val) {
			dofTypes = val;
			return this;
		}
		
		public Builder addTemplates(Consumer<AtomConnectivity.Builder> block) {
			block.accept(atomConnectivityBuilder);
			return this;
		}
		
		public Builder setResPairCache(ResPairCache val) {
			resPairCache = val;
			return this;
		}
		
		public EnergyCalculator build() {
			
			// if no explicit type was picked, pick the best one now
			if (type == null) {
				if (parallelism.numGpus > 0) {
					type = Type.pickBest(dofTypes);
				} else {
					type = Type.Cpu;
				}
			}
			
			// make a res pair cache if needed
			if (resPairCache == null) {
				AtomConnectivity connectivity = atomConnectivityBuilder
					.setParallelism(Parallelism.makeCpu(parallelism.numThreads))
					.build();
				resPairCache = new ResPairCache(ffparams, connectivity);
			}
			
			return new EnergyCalculator(parallelism, type, resPairCache);
		}
	}
	
	public static enum Type {
		
		CpuOriginalCCD {
			
			@Override
			public boolean isSupported() {
				return true;
			}
			
			@Override
			public Context makeContext(Parallelism parallelism, ResPairCache resPairCache) {
				
				return new Context() {{
					EnergyFunctionGenerator egen = new EnergyFunctionGenerator(resPairCache.ffparams);
					numStreams = parallelism.numThreads;
					efuncs = (interactions, mol) -> egen.interactionEnergy(new ForcefieldInteractions(interactions, mol));
					minimizers = (f) -> new CCDMinimizer(f, false);
				}};
			}
		},
		Cpu {
			
			@Override
			public boolean isSupported() {
				return true;
			}
			
			@Override
			public Context makeContext(Parallelism parallelism, ResPairCache resPairCache) {
				
				return new Context() {{
					numStreams = parallelism.numThreads;
					efuncs = (interactions, mol) -> new ResidueForcefieldEnergy(resPairCache, interactions, mol);
					minimizers = (f) -> new SimpleCCDMinimizer(f);
				}};
			}
		},
		Cuda {
			
			@Override
			public boolean isSupported() {
				return !edu.duke.cs.osprey.gpu.cuda.Gpus.get().getGpus().isEmpty();
			}
			
			@Override
			public Context makeContext(Parallelism parallelism, ResPairCache resPairCache) {
				
				// use the Cuda GPU energy function, but do CCD on the CPU
				// (the GPU CCD implementation can't handle non-dihedral dofs yet)
				return new Context() {
					
					private GpuStreamPool pool;
					
					{
						pool = new GpuStreamPool(parallelism.numGpus, parallelism.numStreamsPerGpu);
						numStreams = pool.getNumStreams();
						efuncs = (interactions, mol) -> new GpuForcefieldEnergy(resPairCache.ffparams, new ForcefieldInteractions(interactions, mol), pool);
						minimizers = (f) -> new SimpleCCDMinimizer(f);
						needsCleanup = true;
					}
					
					@Override
					public void cleanup() {
						pool.cleanup();
						needsCleanup = false;
					}
				};
			}
		},
		ResidueCuda { // TODO: eventually replace cuda?
			
			@Override
			public boolean isSupported() {
				return !edu.duke.cs.osprey.gpu.cuda.Gpus.get().getGpus().isEmpty();
			}
			
			@Override
			public Context makeContext(Parallelism parallelism, ResPairCache resPairCache) {
				
				// use the Cuda GPU energy function, but do CCD on the CPU
				// (the GPU CCD implementation can't handle non-dihedral dofs yet)
				return new Context() {
					
					private GpuStreamPool pool;
					
					{
						pool = new GpuStreamPool(parallelism.numGpus, parallelism.numStreamsPerGpu);
						numStreams = pool.getNumStreams();
						efuncs = (interactions, mol) -> new ResidueForcefieldEnergyCuda(pool, resPairCache, interactions, mol);
						minimizers = (f) -> new SimpleCCDMinimizer(f);
						needsCleanup = true;
					}
					
					@Override
					public void cleanup() {
						pool.cleanup();
						needsCleanup = false;
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
			public Context makeContext(Parallelism parallelism, ResPairCache resPairCache) {
				
				// use a CPU energy function, but send it to the Cuda CCD minimizer
				// (which has a built-in GPU energy function)
				return new Context() {
					
					private GpuStreamPool pool;
					
					{
						pool = new GpuStreamPool(parallelism.numGpus, parallelism.numStreamsPerGpu);
						numStreams = pool.getNumStreams();
						efuncs = (interactions, mol) -> new BigForcefieldEnergy(resPairCache.ffparams, new ForcefieldInteractions(interactions, mol), BufferTools.Type.Direct);
						minimizers = (mof) -> new CudaCCDMinimizer(pool, mof);
						needsCleanup = true;
					}
					
					@Override
					public void cleanup() {
						pool.cleanup();
						needsCleanup = false;
					}
				};
			}
		},
		ResidueCudaCCD {
			
			@Override
			public boolean isSupported() {
				return ResidueCuda.isSupported();
			}
			
			@Override
			public Context makeContext(Parallelism parallelism, ResPairCache resPairCache) {
				
				// use a CPU energy function, but send it to the Cuda CCD minimizer
				// (which has a built-in GPU energy function)
				return new Context() {
					
					private GpuStreamPool pool;
					
					{
						pool = new GpuStreamPool(parallelism.numGpus, parallelism.numStreamsPerGpu);
						numStreams = pool.getNumStreams();
						efuncs = (interactions, mol) -> new ResidueForcefieldEnergy(resPairCache, interactions, mol);
						minimizers = (mof) -> new ResidueCudaCCDMinimizer(pool, mof);
						needsCleanup = true;
					}
					
					@Override
					public void cleanup() {
						pool.cleanup();
						needsCleanup = false;
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
			public Context makeContext(Parallelism parallelism, ResPairCache resPairCache) {
				
				// use the CPU CCD minimizer, with an OpenCL energy function
				return new Context() {
					
					private GpuQueuePool pool;
					
					{
						pool = new GpuQueuePool(parallelism.numGpus, parallelism.numStreamsPerGpu);
						numStreams = pool.getNumQueues();
						efuncs = (interactions, mol) -> new GpuForcefieldEnergy(resPairCache.ffparams, new ForcefieldInteractions(interactions, mol), pool);
						minimizers = (mof) -> new SimpleCCDMinimizer(mof);
						needsCleanup = true;
					}
					
					@Override
					public void cleanup() {
						pool.cleanup();
						needsCleanup = false;
					}
				};
			}
		};
		
		public static interface EfuncFactory {
			EnergyFunction make(ResidueInteractions inters, Molecule mol);
		}
		
		public static abstract class Context {
			
			public int numStreams;
			public EfuncFactory efuncs;
			public Factory<Minimizer,ObjectiveFunction> minimizers;
			
			protected boolean needsCleanup = false;
			
			public void cleanup() {
				// do nothing by default
			}
			
			@Override
			protected void finalize()
			throws Throwable {
				try {
					if (needsCleanup) {
						System.err.println("WARNING: " + getClass().getName() + " was garbage collected, but not cleaned up. Attempting cleanup now");
						cleanup();
					}
				} finally {
					super.finalize();
				}
			}
		}
		
		public abstract boolean isSupported();
		public abstract Context makeContext(Parallelism parallelism, ResPairCache resPairCache);
		
		public static Type pickBest(DofTypes dofTypes) {
			
			// only use residue functions automatically (instead of older ones),
			// since they support residue pair weights and offsets and the older ones don't
			if (ResidueCuda.isSupported()) {
				switch (dofTypes) {
					case None:
					case OnlyDihedrals: return ResidueCudaCCD;
					case Any: return ResidueCuda;
				}
			}
			
			// NOTE: OpenCL kernel doesn't support residue pair weights and offsets, so don't pick it automatically
			
			// fallback to CPU, it's always supported
			return Cpu;
		}
	}
	
	public final Parallelism parallelism;
	public final TaskExecutor tasks;
	public final Type type;
	public final Type.Context context;
	public final ResPairCache resPairCache;
	
	private EnergyCalculator(Parallelism parallelism, Type type, ResPairCache resPairCache) {
		
		this.parallelism = parallelism;
		this.tasks = parallelism.makeTaskExecutor();
		this.type = type;
		this.resPairCache = resPairCache;
		
		context = type.makeContext(parallelism, resPairCache);
	}
	
	@Override
	public void clean() {
		context.cleanup();
		tasks.clean();
	}
	
	/**
	 * Calculate the energy of a molecule. If the molecule has continuous degrees of freedom,
	 * they will be minimized within the specified bounds before calculating the energy.
	 * 
	 * @param pmol The molecule
	 * @param bounds Bounds for continuous degrees of freedom for the minimization, if any
	 * @param inters Residue interactions for the energy function
	 * @return The calculated energy
	 */
	public double calcEnergy(ParametricMolecule pmol, DofBounds bounds, ResidueInteractions inters) {
		
		// short circuit: no inters, no energy!
		if (inters.size() <= 0) {
			return 0;
		}
		
		// get the energy function
		EnergyFunction efunc = context.efuncs.make(inters, pmol.mol);
		try {
			
			// get the energy
			double energy;
			if (bounds.size() > 0) {
				
				// minimize it
				Minimizer minimizer = context.minimizers.make(new MoleculeObjectiveFunction(pmol, bounds, efunc));
				try {
					energy = minimizer.minimize().energy;
				} finally {
					Minimizer.Tools.cleanIfNeeded(minimizer);
				}
				
			} else {
				
				// otherwise, just use the score
				energy = efunc.getEnergy();
			}
			
			return energy;
			
		} finally {
			EnergyFunction.Tools.cleanIfNeeded(efunc);
		}
	}
        
        
        public MoleculeObjectiveFunction makeEnergyObjFcn(ParametricMolecule pmol, DofBounds bounds, ResidueInteractions inters) {
                //represent the energy indicated by inters as a function of the degrees of freedom in pmol,
                //valid within the specified bounds
                //calcEnergy gives the minimum of this function
		EnergyFunction efunc = context.efuncs.make(inters, pmol.mol);
		return new MoleculeObjectiveFunction(pmol, bounds, efunc);
	}
        
        
        public void writeMinimizedStruct(ParametricMolecule pmol, DofBounds bounds, ResidueInteractions inters, String fileName) {
		
		// short circuit: no inters, no energy!
		if (inters.size() <= 0) {
			throw new RuntimeException("ERROR: Can't minimize struct with no energy");
		}
		
		// get the energy function
		EnergyFunction efunc = context.efuncs.make(inters, pmol.mol);
		try {
			
			if (bounds.size() > 0) {
				
				// minimize it
				Minimizer minimizer = context.minimizers.make(new MoleculeObjectiveFunction(pmol, bounds, efunc));
				try {
					double energy = minimizer.minimize().energy;
                                        PDBIO.writeFile(pmol.mol, null, energy, fileName);
				} finally {
					Minimizer.Tools.cleanIfNeeded(minimizer);
				}
				
			} else {
                            //no need to minimize
                            double energy = efunc.getEnergy();
                            PDBIO.writeFile(pmol.mol, null, energy, fileName);
			}
		} finally {
			EnergyFunction.Tools.cleanIfNeeded(efunc);
		}
	}
        
        
	
	/**
	 * Asynchronous version of {@link #calcEnergy(ParametricMolecule,DofBouds,ResidueInteractions)}.
	 * 
	 * @param pmol The molecule
	 * @param bounds Bounds for continuous degrees of freedom for the minimization, if any
	 * @param inters Residue interactions for the energy function
	 * @param listener Callback function that will receive the energy. The callback is called on a
	 *                 listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(ParametricMolecule pmol, DofBounds bounds, ResidueInteractions inters, TaskListener<Double> listener) {
		tasks.submit(() -> calcEnergy(pmol, bounds, inters), listener);
	}
}
