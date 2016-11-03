package edu.duke.cs.osprey.energy.forcefield;

import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import com.jogamp.opencl.CLEvent;
import com.jogamp.opencl.CLEvent.CommandType;
import com.jogamp.opencl.CLEvent.ProfilingCommand;
import com.jogamp.opencl.CLException;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.ForcefieldKernel;
import edu.duke.cs.osprey.gpu.cuda.Context;
import edu.duke.cs.osprey.gpu.cuda.ContextPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.ForcefieldKernelOneBlockCuda;
import edu.duke.cs.osprey.gpu.opencl.GpuQueue;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.gpu.opencl.ProfilingEvents;
import edu.duke.cs.osprey.gpu.opencl.kernels.ForcefieldKernelOpenCL;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;

public class GpuForcefieldEnergy implements EnergyFunction.DecomposableByDof, EnergyFunction.NeedsCleanup {
	
	private static final long serialVersionUID = -9142317985561910731L;
	
	private class KernelBuilder {
		
		private ForcefieldKernel kernel;
		
		public KernelBuilder() {
			kernel = null;
		}
		
		public ForcefieldKernel get() {
			
			// do we need to rebuild the forcefield?
			// NOTE: make sure we check for chemical changes even if the kernel is null
			// so we acknowledge any chemical changes that may have happened 
			if (hasChemicalChanges() || kernel == null) {
				
				// cleanup any old kernel if needed
				if (kernel != null) {
					kernel.cleanup();
					kernel = null;
				}
				
				try {
					
					// prep the kernel, upload precomputed data
					if (openclQueue != null) {
						kernel = new ForcefieldKernelOpenCL(openclQueue);
						kernel.setForcefield(new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct));
					} else if (cudaContext != null) {
						// TEMP
						//kernel = new ForcefieldKernelCuda(cudaContext);
						kernel = new ForcefieldKernelOneBlockCuda(cudaContext, new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct));
					} else {
						throw new Error("bad gpu queue/context configuration, this is a bug");
					}
					
				} catch (IOException ex) {
					
					// if we can't find the gpu kernel source, that's something a programmer needs to fix
					throw new Error("can't initialize gpu kernel", ex);
				}
			}
			
			return kernel;
		}
		
		private boolean hasChemicalChanges() {
			
			// look for residue template changes so we can rebuild the forcefield
			boolean hasChanges = false;
			for (AtomGroup[] pair : interactions) {
				for (AtomGroup group : pair) {
					if (group.hasChemicalChange()) {
						hasChanges = true;
					}
					group.ackChemicalChange();
				}
			}
			return hasChanges;
		}
		
		public void cleanup() {
			if (kernel != null) {
				kernel.cleanup();
				kernel = null;
			}
		}
	}
	
	private ForcefieldParams ffparams;
	private ForcefieldInteractions interactions;
	private GpuQueuePool openclQueuePool;
	private GpuQueue openclQueue;
	private ContextPool cudaContextPool;
	private Context cudaContext;
	private KernelBuilder kernelBuilder;
	private BigForcefieldEnergy.Subset ffsubset;
	
	private int numProfilingRuns;
	private Stopwatch uploadStopwatch;
	private Stopwatch kernelStopwatch;
	private Stopwatch downloadStopwatch;
	
	private GpuForcefieldEnergy() {
		numProfilingRuns = 0;
		uploadStopwatch = null;
		kernelStopwatch = null;
		downloadStopwatch = null;
	}
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions, GpuQueuePool queuePool) {
		this();
		this.ffparams = ffparams;
		this.interactions = interactions;
		this.openclQueuePool = queuePool;
		this.openclQueue = queuePool.checkout();
		this.cudaContextPool = null;
		this.cudaContext = null;
		this.kernelBuilder = new KernelBuilder();
		this.ffsubset = null;
	}
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions, ContextPool cudaContextPool) {
		this();
		this.ffparams = ffparams;
		this.interactions = interactions;
		this.openclQueuePool = null;
		this.openclQueue = null;
		this.cudaContextPool = cudaContextPool;
		this.cudaContext = cudaContextPool.checkout();
		this.kernelBuilder = new KernelBuilder();
		this.ffsubset = null;
	}
	
	public GpuForcefieldEnergy(GpuForcefieldEnergy parent, ForcefieldInteractions interactions) {
		this();
		this.ffparams = parent.ffparams;
		this.interactions = interactions;
		if (parent.openclQueue != null) {
			this.openclQueuePool = null;
			this.openclQueue = parent.openclQueue;
			this.cudaContextPool = null;
			this.cudaContext = null;
		} else {
			this.openclQueuePool = null;
			this.openclQueue = null;
			this.cudaContextPool = null;
			this.cudaContext = parent.cudaContext;
		}
		this.kernelBuilder = parent.kernelBuilder;
		this.ffsubset = getKernel().getForcefield().new Subset(interactions);
	}
	
	public ForcefieldKernel getKernel() {
		return kernelBuilder.get();
	}
	
	public BigForcefieldEnergy.Subset getSubset() {
		if (ffsubset != null) {
			return ffsubset;
		}
		return getKernel().getForcefield().getFullSubset();
	}
	
	public void startProfile(int numRuns) {
		numProfilingRuns = numRuns;
		if (getKernel() instanceof ForcefieldKernelOpenCL) {
			((ForcefieldKernelOpenCL)getKernel()).setProfilingEvents(new ProfilingEvents(numRuns*3));
		}
		uploadStopwatch = new Stopwatch();
		kernelStopwatch = new Stopwatch();
		downloadStopwatch = new Stopwatch();
	}
	
	public String dumpProfile() {
		
		StringBuilder buf = new StringBuilder();
		
		if (getKernel() instanceof ForcefieldKernelOpenCL) {
			ForcefieldKernelOpenCL kernel = (ForcefieldKernelOpenCL)getKernel();
			
			// dump gpu profile
			buf.append("GPU Profile:\n");
			Map<CommandType,Long> sums = new EnumMap<>(CommandType.class);
			ProfilingEvents events = kernel.getProfilingEvents();
			for (CLEvent event : events.getCLEvents()) {
				try {
					long startNs = event.getProfilingInfo(ProfilingCommand.START);
					long endNs = event.getProfilingInfo(ProfilingCommand.END);
					long diffNs = endNs - startNs;
					CommandType type = event.getType();
					if (sums.containsKey(type)) {
						sums.put(type, sums.get(type) + diffNs);
					} else {
						sums.put(type, diffNs);
					}
				} catch (CLException.CLProfilingInfoNotAvailableException ex) {
					// don't care, just skip it
				}
			}
			for (Map.Entry<CommandType,Long> entry : sums.entrySet()) {
				buf.append(String.format(
					"\t%16s %s\n",
					entry.getKey(),
					TimeFormatter.format(entry.getValue()/numProfilingRuns, TimeUnit.MICROSECONDS)
				));
			}
			events.cleanup();
			kernel.setProfilingEvents(null);
		}
		
		// dump cpu profile
		buf.append("CPU Profile:\n");
		buf.append(String.format("\t%16s %s\n", "Upload", TimeFormatter.format(uploadStopwatch.getTimeNs()/numProfilingRuns, TimeUnit.MICROSECONDS)));
		buf.append(String.format("\t%16s %s\n", "Kernel", TimeFormatter.format(kernelStopwatch.getTimeNs()/numProfilingRuns, TimeUnit.MICROSECONDS)));
		buf.append(String.format("\t%16s %s\n", "Download", TimeFormatter.format(downloadStopwatch.getTimeNs()/numProfilingRuns, TimeUnit.MICROSECONDS)));
		
		numProfilingRuns = 0;
		
		return buf.toString();
	}
	
	@Override
	public double getEnergy() {
		
		// PROFILING
		//Profiler p = new Profiler("upload");
		
		boolean isProfiling = numProfilingRuns > 0;
		if (isProfiling) {
			uploadStopwatch.resume();
		}
		
		// upload data
		ForcefieldKernel kernel = getKernel();
		kernel.setSubset(getSubset());
		kernel.uploadCoordsAsync();
		
		// PROFILING
		//kernel.waitForGpu();
		//p.start("kernel");
		
		if (isProfiling) {
			kernel.waitForGpu();
			uploadStopwatch.stop();
			kernelStopwatch.resume();
		}
		
		// compute the energies
		kernel.runAsync();
		
		// PROFILING
		//kernel.waitForGpu();
		//p.start("download");
		
		if (isProfiling) {
			kernel.waitForGpu();
			kernelStopwatch.stop();
			downloadStopwatch.resume();
		}
		
		// read the results
		double energy = kernel.downloadEnergySync();
		
		// PROFILING
		//System.out.println(p.makeReport(TimeUnit.MICROSECONDS));
		
		if (isProfiling) {
			downloadStopwatch.stop();
		}
		
		return energy;
	}
	
	@Override
	public void cleanup() {
		
		kernelBuilder.cleanup();
		
		if (openclQueuePool != null) {
			openclQueuePool.release(openclQueue);
		}
		
		if (cudaContextPool != null) {
			cudaContextPool.release(cudaContext);
		}
	}

	@Override
	public List<EnergyFunction> decomposeByDof(Molecule m, List<DegreeOfFreedom> dofs) {
		
		List<EnergyFunction> efuncs = new ArrayList<>();
		Map<Residue,GpuForcefieldEnergy> efuncCache = new HashMap<>();
		
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
					efunc = new GpuForcefieldEnergy(this, interactions.makeSubsetByResidue(res));
					efuncCache.put(res, efunc);
				}
				efuncs.add(efunc);
			}
		}
		
		return efuncs;
	}
}
