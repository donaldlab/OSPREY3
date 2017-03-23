package edu.duke.cs.osprey.minimization;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Factory;

public class ExampleParallelMinimization {
	
	@SuppressWarnings("unused")
	private static List<EnergiedConf> minimize(ForcefieldParams ffparams, ConfSpace confSpace, List<Residue> shellResidues, List<ScoredConf> confs) {
		
		// get forcefield interactions
		Factory<ForcefieldInteractions,Molecule> ffinteractions;
		
		// use full conf interactions
		ffinteractions = (mol) -> FFInterGen.makeFullConf(confSpace, shellResidues, mol);
		
		// or full molecule interactions
		ffinteractions = (mol) -> FFInterGen.makeFullMol(mol);
		
		// build the minimizer
		ConfMinimizer minimizer;
		
		// use the CPU
		minimizer = new CpuConfMinimizer.Builder(ffparams, ffinteractions, confSpace).build();
		
		// or use multiple CPU threads
		minimizer = new CpuConfMinimizer.Builder(ffparams, ffinteractions, confSpace).setNumThreads(4).build();
		
		// or use the GPU (will automatically pick  best implementation)
		minimizer = new GpuConfMinimizer.Builder(ffparams, ffinteractions, confSpace).build();
		
		// or use multiple GPUs, or multiple streams
		minimizer = new GpuConfMinimizer.Builder(ffparams, ffinteractions, confSpace).setGpuInfo(GpuConfMinimizer.Type.OpenCL, 8, 2).build();
		
		// or explicitly pick the super duper fast Cuda CCD minimizer and give it lots of streams
		minimizer = new GpuConfMinimizer.Builder(ffparams, ffinteractions, confSpace).setGpuInfo(GpuConfMinimizer.Type.CudaCCD, 8, 16).build();
		
		// do the minimization (in parallel, if supported)
		List<EnergiedConf> minimizedConfs = minimizer.minimize(confs);
		
		// don't forget to clean up when you're done
		// (to release gpu resources, stop worker threads, etc)
		minimizer.cleanup();
		
		return minimizedConfs;
	}
}
