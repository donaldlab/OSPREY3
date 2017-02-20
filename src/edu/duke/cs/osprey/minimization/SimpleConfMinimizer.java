package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.control.Defaults;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gmec.ConfEnergyCalculator;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.Progress;

/**
 * Computes the energy of a conformation using the desired forcefield.
 * 
 * If a conformation has continuous degrees of freedom, minimization will be performed before forcefield evaluation.
 */
public class SimpleConfMinimizer implements ConfEnergyCalculator.Async {
	
	public static class Builder {
		
		private SimpleConfSpace confSpace;
		private Parallelism parallelism;
		private Type type;
		
		/**
		 * Is the number of conformations to be minimized unknown in advance?
		 * 
		 * @todo describe conf streaming and ThreadPoolTaskExecutor buffers. or just get rid of the buffering entirely.
		 */
		private boolean isStreaming = false;
		private ForcefieldParams ffparams = Defaults.forcefieldParams;
		
		public Builder(SimpleConfSpace confSpace) {
			this.confSpace = confSpace;
			this.parallelism = Parallelism.makeDefault();
			this.type = null;
		}
		
		public Builder setParallelism(Parallelism val) {
			parallelism = val;
			return this;
		}
		
		public Builder setType(Type val) {
			type = val;
			return this;
		}
		
		public Builder setStreaming(boolean val) {
			isStreaming = val;
			return this;
		}
		
		public Builder setForcefieldParams(ForcefieldParams val) {
			ffparams = val;
			return this;
		}
		
		public SimpleConfMinimizer build() {
			
			// if no explict type was picked, pick the best one now
			if (type == null) {
				if (parallelism.numGpus > 0) {
					type = Type.pickBest(confSpace);
				} else {
					type = Type.Cpu;
				}
			}
			
			return new SimpleConfMinimizer(
				confSpace,
				parallelism,
				isStreaming,
				type,
				ffparams
			);
		}
	}
	
	public static Builder builder(SimpleConfSpace confSpace) {
		return new Builder(confSpace);
	}
	
	public static enum Type {
		
		Cpu {
			
			@Override
			public boolean isSupported() {
				return true;
			}
			
			@Override
			public Context makeContext(Parallelism parallelism, ForcefieldParams ffparams) {
				
				return new Context() {
					
					private EnergyFunctionGenerator egen;
					
					{
						egen = new EnergyFunctionGenerator(ffparams);
						numStreams = parallelism.numThreads;
						efuncs = (interactions) -> egen.interactionEnergy(interactions);
						minimizers = (f) -> new SimpleCCDMinimizer(f);
					}
				};
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
	public final boolean isStreaming;
	public final TaskExecutor tasks;
	public final Type type;
	public final ForcefieldParams ffparams;
	
	private Type.Context context;
	
	private SimpleConfMinimizer(SimpleConfSpace confSpace, Parallelism parallelism, boolean isStreaming, Type type, ForcefieldParams ffparams) {
		
		this.confSpace = confSpace;
		this.parallelism = parallelism;
		this.isStreaming = isStreaming;
		this.tasks = parallelism.makeTaskExecutor(isStreaming);
		this.type = type;
		this.ffparams = ffparams;
		
		context = type.makeContext(parallelism, ffparams);
	}
	
	public void cleanup() {
		parallelism.cleanupTaskExecutor(tasks);
		context.cleanup();
	}
	
	public EnergiedConf minimizeSync(ScoredConf conf) {
		return new EnergiedConf(conf, minimizeSync(conf.getAssignments()).energy);
	}
	
	public Minimizer.Result minimizeSync(int[] conf) {
		return minimizeSync(new RCTuple(conf));
	}
	
	public Minimizer.Result minimizeSync(RCTuple conf) {
		
		// make the mol in the conf
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		
		EnergyFunction efunc = null;
		Minimizer minimizer = null;
		
		try {
			
			// get the energy function
			efunc = context.efuncs.make(FFInterGen.makeFullConf(confSpace, pmol.mol));
			
			DofBounds bounds = confSpace.makeBounds(conf);
			if (bounds.size() > 0) {
				
				// minimize it
				minimizer = context.minimizers.make(new MoleculeObjectiveFunction(pmol, bounds, efunc));
				return minimizer.minimize();
				
			} else {
				
				// otherwise, just use the score
				return new Minimizer.Result(null, efunc.getEnergy());
			}
			
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
	
	public static interface Listener<T> {
		void onMinimized(T result);
	}
	
	public void minimizeAsync(ScoredConf conf) {
		minimizeAsync(conf, null);
	}
	
	public void minimizeAsync(ScoredConf conf, Listener<EnergiedConf> listener) {
		
		class Task implements Runnable {
			
			public EnergiedConf econf;
			
			@Override
			public void run() {
				econf = minimizeSync(conf);
			}
		}
		
		// submit the minimization task and chain the listener if needed
		if (listener == null) {
			tasks.submit(new Task());
		} else {
			tasks.submit(new Task(), (task) -> {
				listener.onMinimized(task.econf);
			});
		}
	}
	
	public void minimizeAsync(int[] conf) {
		minimizeAsync(conf, null);
	}
	
	public void minimizeAsync(int[] conf, Listener<Minimizer.Result> listener) {
		minimizeAsync(new RCTuple(conf), listener);
	}
	
	public void minimizeAsync(RCTuple conf) {
		minimizeAsync(conf, null);
	}
	
	public void minimizeAsync(RCTuple conf, Listener<Minimizer.Result> listener) {
		
		class Task implements Runnable {
			
			public Minimizer.Result result;
			
			@Override
			public void run() {
				result = minimizeSync(conf);
			}
		}
		
		// submit the minimization task and chain the listener if needed
		if (listener == null) {
			tasks.submit(new Task());
		} else {
			tasks.submit(new Task(), (task) -> {
				listener.onMinimized(task.result);
			});
		}
	}
	
	@Override
	public int getParallelism() {
		return parallelism.getParallelism();
	}
	
	@Override
	public void waitForSpace() {
		tasks.waitForSpace();
	}
	
	@Override
	public void waitForFinish() {
		tasks.waitForFinish();
	}
	
	@Override
	public EnergiedConf calcEnergy(ScoredConf conf) {
		return minimizeSync(conf);
	}

	@Override
	public void calcEnergyAsync(ScoredConf conf, ConfEnergyCalculator.Async.Listener listener) {
		minimizeAsync(conf, (econf) -> {
			listener.onEnergy(econf);
		});
	}
	
	public List<EnergiedConf> minimizeSync(List<ScoredConf> confs) {
		return minimizeSync(confs, false);
	}
	
	public List<EnergiedConf> minimizeSync(List<ScoredConf> confs, boolean reportProgress) {
		
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
			
			minimizeAsync(confs.get(i), (econf) -> {
				
				// save the minimized energy
				econfs.set(fi, econf);
				
				// update progress if needed
				if (progress != null) {
					progress.incrementProgress();
				}
			});
		}
		waitForFinish();
		
		return econfs;
	}
}
