package edu.duke.cs.osprey.energy;

import java.io.File;
import java.util.function.Consumer;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
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
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
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

		/** True to minimize continuous degrees of freedom in conformations. False to use only rigid structures. */
		private boolean isMinimizing = true;

		/**
		 * When set, if a minimized energy is below this threshold, assume the conformation
		 * fell into an infinite well due to a severe clash (this is a known flaw of the AMBER
		 * forcefield, all parameterizations). Then attempt the minimization again,
		 * but try to resolve the clash first with a purely vdW forcefield. If the clash cannot
		 * be resolved, and the conformation falls into the well again, return positive infinity
		 * energy instead of effectively negative infinity.
		 *
		 * If you wish to use the inifinite well detection and resolution, -10,000 kcal/mol is
		 * probably a good threshold to choose.
		 */
		private Double infiniteWellEnergy = null;

		/**
		 * If set, the minimizer will always use a simpler vdW forcefield to attempt to resolve clashes
		 * before applying the full forcefield. If the energy does not minimize to below the threshold,
		 * the full forcefield will not be used, and a positive infinity energy will be returned instead.
		 */
		private Double alwaysResolveClashesEnergy = null;

		public Builder(ResidueTemplateLibrary templateLib, ForcefieldParams ffparams) {
			this.ffparams = ffparams;
			this.atomConnectivityBuilder.addTemplates(templateLib);
		}

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

		public Builder setIsMinimizing(boolean val) {
			this.isMinimizing = val;
			return this;
		}

		public Builder setInfiniteWellEnergy(Double val) {
			infiniteWellEnergy = val;
			return this;
		}

		public Builder setAlwaysResolveClashesEnergy(Double val) {
			alwaysResolveClashesEnergy = val;
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
			
			return new EnergyCalculator(parallelism, type, resPairCache, isMinimizing, infiniteWellEnergy, alwaysResolveClashesEnergy);
		}
	}

	/**
	 * Create an energy calculator that shares resources with an existing energy calculator,
	 * (including thread pools, GPUs, atom connectivity cache, forcefield),
	 * but can have different configuration.
	 */
	public static class SharedBuilder {

		private EnergyCalculator parent;

		private boolean isMinimizing = true;

		public SharedBuilder(EnergyCalculator parent) {
			this.parent = parent;
		}

		public SharedBuilder setIsMinimizing(boolean val) {
			this.isMinimizing = val;
			return this;
		}

		public EnergyCalculator build() {
			return new EnergyCalculator(parent, isMinimizing) {

				@Override
				public void clean() {
					// ignore cleanup calls, since we're sharing resources with our parent
				}
			};
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

	public static class EnergiedParametricMolecule {

		public final ParametricMolecule pmol;
		public final ResidueInteractions inters;
		public final DoubleMatrix1D params;
		public final double energy;

		public EnergiedParametricMolecule(ParametricMolecule pmol, ResidueInteractions inters, double energy) {
			this(pmol, inters, null, energy);
		}

		public EnergiedParametricMolecule(ParametricMolecule pmol, ResidueInteractions inters, DoubleMatrix1D params, double energy) {
			this.pmol = pmol;
			this.inters = inters;
			this.params = params;
			this.energy = energy;
		}
	}
	
	public final Parallelism parallelism;
	public final TaskExecutor tasks;
	public final Type type;
	public final Type.Context context;
	public final ResPairCache resPairCache;
	public final boolean isMinimizing;
	public final Double infiniteWellEnergy;
	public final Double alwaysResolveClashesEnergy;

	private final Type.Context cpuContext; // for vdW forcefields
	
	private EnergyCalculator(Parallelism parallelism, Type type, ResPairCache resPairCache, boolean isMinimizing, Double infiniteWellEnergy, Double alwaysResolveClashesEnergy) {

		this.parallelism = parallelism;
		this.tasks = parallelism.makeTaskExecutor();
		this.type = type;
		this.context = type.makeContext(parallelism, resPairCache);
		this.resPairCache = resPairCache;
		this.isMinimizing = isMinimizing;
		this.infiniteWellEnergy = infiniteWellEnergy;
		this.alwaysResolveClashesEnergy = alwaysResolveClashesEnergy;

		// make a CPU context if we need to do vdW forcefields
		// TODO: implement vdW forcefield on the GPU too?
		if (infiniteWellEnergy != null || alwaysResolveClashesEnergy != null) {
			this.cpuContext = Type.Cpu.makeContext(parallelism, resPairCache);
		} else {
			this.cpuContext = null;
		}
	}

	private EnergyCalculator(EnergyCalculator parent, boolean isMinimizing) {

		this.parallelism = parent.parallelism;
		this.tasks = parent.tasks;
		this.type = parent.type;
		this.context = parent.context;
		this.resPairCache = parent.resPairCache;
		this.isMinimizing = isMinimizing;
		this.infiniteWellEnergy = parent.infiniteWellEnergy;
		this.alwaysResolveClashesEnergy = parent.alwaysResolveClashesEnergy;

		this.cpuContext = parent.cpuContext;
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
	 * @param inters Residue interactions for the energy function
	 * @return The calculated energy and the associated molecule pose
	 */
	public EnergiedParametricMolecule calcEnergy(ParametricMolecule pmol, ResidueInteractions inters) {
		
		// short circuit: no inters, no energy!
		if (inters.size() <= 0) {
			return new EnergiedParametricMolecule(pmol, null, 0);
		}

		// if we don't need to minimize, just return the energy of the current pose
		if (!isMinimizing || pmol.dofBounds.size() <= 0) {
			try (EnergyFunction efunc = context.efuncs.make(inters, pmol.mol)) {
				return new EnergiedParametricMolecule(pmol, inters, null, efunc.getEnergy());
			}
		}

		// start at the center of the voxel
		DoubleMatrix1D x = DoubleFactory1D.dense.make(pmol.dofs.size());
		pmol.dofBounds.getCenter(x);

		if (alwaysResolveClashesEnergy != null) {

			Minimizer.Result vdwResult = minimizeWithVdw(pmol, inters, x);

			// if we didn't resolve the clash, return +inf energy
			if (vdwResult.energy >= alwaysResolveClashesEnergy) {
				return new EnergiedParametricMolecule(pmol, inters, vdwResult.dofValues, Double.POSITIVE_INFINITY);
			}

			x = vdwResult.dofValues;
		}

		// minimize using the full forcefield
		try (EnergyFunction efunc = context.efuncs.make(inters, pmol.mol)) {
			MoleculeObjectiveFunction mof = new MoleculeObjectiveFunction(pmol, efunc);
			try (Minimizer minimizer = context.minimizers.make(mof)) {

				Minimizer.Result result = minimizer.minimizeFrom(x);

				// did we fall into an infinite energy well?
				if (isInfiniteWell(result.energy)) {

					// try to resolve the clash and try the minimization again
					Minimizer.Result vdwResult = minimizeWithVdw(pmol, inters, x);
					result = minimizer.minimizeFrom(vdwResult.dofValues);

					// are we still in the well?
					if (isInfiniteWell(result.energy)) {

						// return positive infinity energy
						return new EnergiedParametricMolecule(pmol, inters, result.dofValues, Double.POSITIVE_INFINITY);
					}

					// we got out of the well, yay!
				}

				return new EnergiedParametricMolecule(pmol, inters, result.dofValues, result.energy);
			}
		}
	}

	private Minimizer.Result minimizeWithVdw(ParametricMolecule pmol, ResidueInteractions inters, DoubleMatrix1D x) {
		try (EnergyFunction efunc = cpuContext.efuncs.make(inters, pmol.mol)) {
			ResidueForcefieldEnergy.Vdw vdwEfunc = new ResidueForcefieldEnergy.Vdw((ResidueForcefieldEnergy)efunc);
			try (Minimizer minimizer = cpuContext.minimizers.make(new MoleculeObjectiveFunction(pmol, vdwEfunc))) {
				return minimizer.minimizeFrom(x);
			}
		}
	}

	private boolean isInfiniteWell(double energy) {
		return infiniteWellEnergy != null && energy <= infiniteWellEnergy;
	}

	public EnergyFunction makeEnergyFunction(EnergiedParametricMolecule epmol) {
		return makeEnergyFunction(epmol.pmol, epmol.inters);
	}

	public EnergyFunction makeEnergyFunction(ParametricMolecule pmol, ResidueInteractions inters) {
		return context.efuncs.make(inters, pmol.mol);
	}

	public MoleculeObjectiveFunction makeEnergyObjFcn(ParametricMolecule pmol, ResidueInteractions inters) {
		//represent the energy indicated by inters as a function of the degrees of freedom in pmol,
		//valid within the specified bounds
		//calcEnergy gives the minimum of this function
		return new MoleculeObjectiveFunction(pmol, makeEnergyFunction(pmol, inters));
	}

	/**
	 * Asynchronous version of {@link #calcEnergy(ParametricMolecule,ResidueInteractions)}.
	 * 
	 * @param pmol The molecule
	 * @param inters Residue interactions for the energy function
	 * @param listener Callback function that will receive the energy and the associated molecule pose.
	 *                 The callback is called on a listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(ParametricMolecule pmol, ResidueInteractions inters, TaskListener<EnergiedParametricMolecule> listener) {
		tasks.submit(() -> calcEnergy(pmol, inters), listener);
	}
}
