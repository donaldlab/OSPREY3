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
import com.jogamp.opencl.CLEventList;
import com.jogamp.opencl.CLException;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.gpu.GpuQueue;
import edu.duke.cs.osprey.gpu.GpuQueuePool;
import edu.duke.cs.osprey.gpu.kernels.ForceFieldKernel;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;

public class GpuForcefieldEnergy implements EnergyFunction.DecomposableByDof, EnergyFunction.NeedsCleanup {
	
	private static final long serialVersionUID = -9142317985561910731L;
	
	private class KernelBuilder {
		
		private ForceFieldKernel.Bound kernel;
		
		public KernelBuilder() {
			kernel = null;
		}
		
		public ForceFieldKernel.Bound get() {
			
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
					kernel = new ForceFieldKernel(queue.getGpu()).bind(queue);
					kernel.setForcefield(new BigForcefieldEnergy(ffparams, interactions, true));
					kernel.uploadStaticAsync();
					
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
	private GpuQueuePool queuePool;
	private GpuQueue queue;
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
		this.queuePool = queuePool;
		this.queue = queuePool.checkout();
		this.kernelBuilder = new KernelBuilder();
		this.ffsubset = null;
	}
	
	public GpuForcefieldEnergy(GpuForcefieldEnergy parent, ForcefieldInteractions interactions) {
		this();
		this.ffparams = parent.ffparams;
		this.interactions = interactions;
		this.queuePool = null;
		this.queue = parent.queue;
		this.kernelBuilder = parent.kernelBuilder;
		this.ffsubset = getKernel().getForcefield().new Subset(interactions);
	}
	
	public boolean isParent() {
		return queuePool != null;
	}
	
	public ForceFieldKernel.Bound getKernel() {
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
		getKernel().initProfilingEvents(numRuns*3);
		uploadStopwatch = new Stopwatch();
		kernelStopwatch = new Stopwatch();
		downloadStopwatch = new Stopwatch();
	}
	
	public String dumpProfile() {
		
		StringBuilder buf = new StringBuilder();
		
		// dump gpu profile
		buf.append("GPU Profile:\n");
		Map<CommandType,Long> sums = new EnumMap<>(CommandType.class);
		ForceFieldKernel.Bound kernel = getKernel();
		CLEventList events = kernel.getProfilingEvents();
		for (CLEvent event : events) {
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
		kernel.clearProfilingEvents();
		
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
		
		boolean isProfiling = numProfilingRuns > 0;
		if (isProfiling) {
			uploadStopwatch.resume();
		}
		
		// upload data
		ForceFieldKernel.Bound kernel = getKernel();
		kernel.setSubset(getSubset());
		kernel.uploadCoordsAsync();
		
		if (isProfiling) {
			kernel.waitForGpu();
			uploadStopwatch.stop();
			kernelStopwatch.resume();
		}
		
		// compute the energies
		kernel.runAsync();
		
		if (isProfiling) {
			kernel.waitForGpu();
			kernelStopwatch.stop();
			downloadStopwatch.resume();
		}
		
		// read the results
		double energy = kernel.downloadEnergySync();
		
		if (isProfiling) {
			downloadStopwatch.stop();
		}
		
		return energy;
	}
	
	@Override
	public void cleanup() {
		if (isParent()) {
			queuePool.release(queue);
			kernelBuilder.cleanup();
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
