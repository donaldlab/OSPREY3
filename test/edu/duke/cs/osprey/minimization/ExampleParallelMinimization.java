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
