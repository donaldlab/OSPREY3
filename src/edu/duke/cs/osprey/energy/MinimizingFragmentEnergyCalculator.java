package edu.duke.cs.osprey.energy;

import java.io.IOException;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
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
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.Profiler;

/**
 * Computes the energy of a molecule fragment using the desired forcefield.
 * 
 * If a fragment has continuous degrees of freedom, minimization will be performed before forcefield evaluation.
 */
public class MinimizingFragmentEnergyCalculator implements FragmentEnergyCalculator.Async {
	
	public static class Builder {
		
		private SimpleConfSpace confSpace;
		private ForcefieldParams ffparams;
		private Parallelism parallelism = Parallelism.makeCpu(1);
		private Type type = null;
		
		public Builder(SimpleConfSpace confSpace, ForcefieldParams ffparams) {
			this.confSpace = confSpace;
			this.ffparams = ffparams;
		}
		
		public Builder setParallelism(Parallelism val) {
			parallelism = val;
			return this;
		}
		
		public Builder setType(Type val) {
			type = val;
			return this;
		}
		
		public MinimizingFragmentEnergyCalculator build() {
			
			// if no explict type was picked, pick the best one now
			if (type == null) {
				if (parallelism.numGpus > 0) {
					type = Type.pickBest(confSpace);
				} else {
					type = Type.Cpu;
				}
			}
			
			return new MinimizingFragmentEnergyCalculator(
				confSpace,
				parallelism,
				type,
				ffparams
			);
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
						efuncs = (interactions, mol) -> {
							try {
								return new ResidueForcefieldEnergyCuda(pool, resPairCache, interactions, mol);
							} catch (IOException ex) {
								throw new Error(ex);
							}
						};
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
		ResidueCudaCCD { // TODO: eventually replace CudaCCD?
			
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
						efuncs = (interactions, mol) -> new ResidueForcefieldEnergy(resPairCache, interactions, mol);
						minimizers = (mof) -> {
							try {
								return new ResidueCudaCCDMinimizer(pool, mof);
							} catch (IOException ex) {
								throw new Error(ex);
							}
						};
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
		
		private static interface EfuncFactory {
			EnergyFunction make(ResidueInteractions inters, Molecule mol);
		}
		
		private static abstract class Context {
			
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
		
		public static Type pickBest(SimpleConfSpace confSpace) {
			
			// prefer cuda over opencl, when both are available
			// only because our cuda code is much better than the opencl code right now
			if (Cuda.isSupported()) {
				if (confSpace.isGpuCcdSupported()) {
					return CudaCCD;
				} else {
					return Cuda;
				}
			}
			
			if (OpenCL.isSupported()) {
				return OpenCL;
			}
			
			// fallback to CPU, it's always supported
			return Cpu;
		}
	}
	
	public final SimpleConfSpace confSpace;
	public final Parallelism parallelism;
	public final TaskExecutor tasks;
	public final Type type;
	
	private Type.Context context;
	
	private MinimizingFragmentEnergyCalculator(SimpleConfSpace confSpace, Parallelism parallelism, Type type, ForcefieldParams ffparams) {
		
		this.confSpace = confSpace;
		this.parallelism = parallelism;
		this.tasks = parallelism.makeTaskExecutor();
		this.type = type;
		
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.setConfSpace(confSpace)
			.setParallelism(parallelism)
			.build();
		ResPairCache resPairCache = new ResPairCache(ffparams, connectivity);
		
		context = type.makeContext(parallelism, resPairCache);
	}
	
	@Override
	public SimpleConfSpace getConfSpace() {
		return confSpace;
	}
	
	public void cleanup() {
		context.cleanup();
	}
	
	@Override
	public double calcEnergy(RCTuple frag, ResidueInteractions inters) {
		
		// TEMP
		Profiler p = new Profiler();
		p.start("mol");
		
		// make the mol in the conf
		ParametricMolecule pmol = confSpace.makeMolecule(frag);
		
		EnergyFunction efunc = null;
		Minimizer minimizer = null;
		
		try {
			
			// TEMP
			p.start("efunc");
			
			// get the energy function
			efunc = context.efuncs.make(inters, pmol.mol);
			
			// TEMP
			p.start("bounds");
				
			// get the energy
			double energy;
			DofBounds bounds = confSpace.makeBounds(frag);
			if (bounds.size() > 0) {
				
				// TEMP
				p.start("mof");
				
				MoleculeObjectiveFunction mof = new MoleculeObjectiveFunction(pmol, bounds, efunc);
				
				// minimize it
				//minimizer = context.minimizers.make(new MoleculeObjectiveFunction(pmol, bounds, efunc));
				
				// TEMP
				p.start("make");
				minimizer = context.minimizers.make(mof);
				
				// TEMP
				p.start("min");
				
				energy = minimizer.minimize().energy;
				
			} else {
				
				// otherwise, just use the score
				energy = efunc.getEnergy();
			}
			
			// TEMP
			p.stop();
			
			return energy;
			
		// make sure we always cleanup the energy function and minimizer
		} finally {
			
			// TEMP
			p.start("cleanup");
			
			if (efunc != null && efunc instanceof EnergyFunction.NeedsCleanup) {
				((EnergyFunction.NeedsCleanup)efunc).cleanup();
			}
			if (minimizer != null && minimizer instanceof Minimizer.NeedsCleanup) {
				((Minimizer.NeedsCleanup)minimizer).cleanup();
			}
			
			// TEMP
			p.stop();
			//System.out.println(p.makeReport(TimeUnit.MILLISECONDS));
		}
	}
	
	@Override
	public void calcEnergyAsync(RCTuple frag, ResidueInteractions inters, FragmentEnergyCalculator.Async.Listener listener) {
		tasks.submit(() -> calcEnergy(frag, inters), listener);
	}
	
	@Override
	public TaskExecutor getTasks() {
		return tasks;
	}
}
