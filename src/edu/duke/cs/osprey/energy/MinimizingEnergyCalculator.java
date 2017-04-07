package edu.duke.cs.osprey.energy;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
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
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.Progress;

/**
 * Computes the energy of a conformation using the desired forcefield.
 * 
 * If a conformation has continuous degrees of freedom, minimization will be performed before forcefield evaluation.
 */
public class MinimizingEnergyCalculator implements ConfEnergyCalculator.Async, FragmentEnergyCalculator.Async {
	
	public static class Builder {
		
		private SimpleConfSpace confSpace;
		private ForcefieldParams ffparams;
		private Parallelism parallelism = Parallelism.makeCpu(1);
		private Type type = null;
		private SimpleReferenceEnergies eref = null;
		
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
		
		public Builder setReferenceEnergies(SimpleReferenceEnergies val) {
			eref = val;
			return this;
		}
		
		public MinimizingEnergyCalculator build() {
			
			// if no explict type was picked, pick the best one now
			if (type == null) {
				if (parallelism.numGpus > 0) {
					type = Type.pickBest(confSpace);
				} else {
					type = Type.Cpu;
				}
			}
			
			return new MinimizingEnergyCalculator(
				confSpace,
				parallelism,
				type,
				ffparams,
				eref
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
			public Context makeContext(Parallelism parallelism, ForcefieldParams ffparams) {
				
				return new Context() {{
					EnergyFunctionGenerator egen = new EnergyFunctionGenerator(ffparams);
					numStreams = parallelism.numThreads;
					efuncs = (interactions) -> egen.interactionEnergy(interactions);
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
			public Context makeContext(Parallelism parallelism, ForcefieldParams ffparams) {
				
				return new Context() {{
					EnergyFunctionGenerator egen = new EnergyFunctionGenerator(ffparams);
					numStreams = parallelism.numThreads;
					efuncs = (interactions) -> egen.interactionEnergy(interactions);
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
			public Context makeContext(Parallelism parallelism, ForcefieldParams ffparams) {
				
				// use the Cuda GPU energy function, but do CCD on the CPU
				// (the GPU CCD implementation can't handle non-dihedral dofs yet)
				return new Context() {
					
					private GpuStreamPool pool;
					
					{
						pool = new GpuStreamPool(parallelism.numGpus, parallelism.numStreamsPerGpu);
						numStreams = pool.getNumStreams();
						efuncs = (interactions) -> new GpuForcefieldEnergy(ffparams, interactions, pool);
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
			public Context makeContext(Parallelism parallelism, ForcefieldParams ffparams) {
				
				// use a CPU energy function, but send it to the Cuda CCD minimizer
				// (which has a built-in GPU energy function)
				return new Context() {
					
					private GpuStreamPool pool;
					
					{
						pool = new GpuStreamPool(parallelism.numGpus, parallelism.numStreamsPerGpu);
						numStreams = pool.getNumStreams();
						efuncs = (interactions) -> new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
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
		OpenCL {
			
			@Override
			public boolean isSupported() {
				return !edu.duke.cs.osprey.gpu.opencl.Gpus.get().getGpus().isEmpty();
			}

			@Override
			public Context makeContext(Parallelism parallelism, ForcefieldParams ffparams) {
				
				// use the CPU CCD minimizer, with an OpenCL energy function
				return new Context() {
					
					private GpuQueuePool pool;
					
					{
						pool = new GpuQueuePool(parallelism.numGpus, parallelism.numStreamsPerGpu);
						numStreams = pool.getNumQueues();
						efuncs = (interactions) -> new GpuForcefieldEnergy(ffparams, interactions, pool);
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
		
		private static abstract class Context {
			
			public int numStreams;
			public Factory<EnergyFunction,ForcefieldInteractions> efuncs;
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
		public abstract Context makeContext(Parallelism parallelism, ForcefieldParams ffparams);
		
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
	public final ForcefieldParams ffparams;
	public final SimpleReferenceEnergies eref;
	
	private Type.Context context;
	
	private MinimizingEnergyCalculator(SimpleConfSpace confSpace, Parallelism parallelism, Type type, ForcefieldParams ffparams, SimpleReferenceEnergies eref) {
		
		this.confSpace = confSpace;
		this.parallelism = parallelism;
		this.tasks = parallelism.makeTaskExecutor();
		this.type = type;
		this.ffparams = ffparams;
		this.eref = eref;
		
		context = type.makeContext(parallelism, ffparams);
	}
	
	public void cleanup() {
		context.cleanup();
	}
	
	@Override
	public double calcEnergy(RCTuple frag, InteractionsFactory intersFactory) {
		
		// make the mol in the conf
		ParametricMolecule pmol = confSpace.makeMolecule(frag);
		
		EnergyFunction efunc = null;
		Minimizer minimizer = null;
		
		try {
			
			// get the energy function
			efunc = context.efuncs.make(intersFactory.make(pmol.mol));
			
			// get the energy
			double energy;
			DofBounds bounds = confSpace.makeBounds(frag);
			if (bounds.size() > 0) {
				
				// minimize it
				minimizer = context.minimizers.make(new MoleculeObjectiveFunction(pmol, bounds, efunc));
				energy = minimizer.minimize().energy;
				
			} else {
				
				// otherwise, just use the score
				energy = efunc.getEnergy();
			}
			
			// apply reference energies if needed
			if (eref != null) {
				energy += eref.getFragmentEnergy(confSpace, frag);
			}
			
			// TODO: entropies
			
			return energy;
			
		// make sure we always cleanup the energy function and minimizer
		} finally {
			if (efunc != null && efunc instanceof EnergyFunction.NeedsCleanup) {
				((EnergyFunction.NeedsCleanup)efunc).cleanup();
			}
			if (minimizer != null && minimizer instanceof Minimizer.NeedsCleanup) {
				((Minimizer.NeedsCleanup)minimizer).cleanup();
			}
		}
	}
	
	@Override
	public void calcEnergyAsync(RCTuple frag, InteractionsFactory intersFactory, FragmentEnergyCalculator.Async.Listener listener) {
		tasks.submit(() -> calcEnergy(frag, intersFactory), listener);
	}
	
	@Override
	public EnergiedConf calcEnergy(ScoredConf conf) {
		RCTuple tup = new RCTuple(conf.getAssignments());
		double energy = calcEnergy(tup, (Molecule mol) -> FFInterGen.makeFullConf(confSpace, mol));
		return new EnergiedConf(conf, energy);
	}

	@Override
	public void calcEnergyAsync(ScoredConf conf, ConfEnergyCalculator.Async.Listener listener) {
		tasks.submit(() -> calcEnergy(conf), listener);
	}
	
	@Override
	public TaskExecutor getTasks() {
		return tasks;
	}
	
	public List<EnergiedConf> calcAllEnergies(List<ScoredConf> confs) {
		return calcAllEnergies(confs, false);
	}
	
	public List<EnergiedConf> calcAllEnergies(List<ScoredConf> confs, boolean reportProgress) {
		
		// allocate space to hold the minimized values
		List<EnergiedConf> econfs = new ArrayList<>(confs.size());
		for (int i=0; i<confs.size(); i++) {
			econfs.add(null);
		}
		
		// track progress if desired
		final Progress progress;
		if (reportProgress) {
			progress = new Progress(confs.size());
		} else {
			progress = null;
		}
		
		// minimize them all
		for (int i=0; i<confs.size(); i++) {
			
			// capture i for the closure below
			final int fi = i;
			
			calcEnergyAsync(confs.get(i), (econf) -> {
				
				// save the minimized energy
				econfs.set(fi, econf);
				
				// update progress if needed
				if (progress != null) {
					progress.incrementProgress();
				}
			});
		}
		tasks.waitForFinish();
		
		return econfs;
	}
}
