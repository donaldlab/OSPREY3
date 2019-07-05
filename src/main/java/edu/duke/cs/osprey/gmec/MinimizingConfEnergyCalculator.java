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

package edu.duke.cs.osprey.gmec;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.ReferenceEnergies;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.minimization.ConfMinimizer;
import edu.duke.cs.osprey.minimization.CpuConfMinimizer;
import edu.duke.cs.osprey.minimization.GpuConfMinimizer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

/** Use the new MinimizingEnergyCalculator instead. Of course, that means you need to switch to the new SimpleConfSpace too. */
@Deprecated
public class MinimizingConfEnergyCalculator implements GMECConfEnergyCalculator.Async {
	
	public static MinimizingConfEnergyCalculator make(ForcefieldParams ffparams, SearchProblem search) {
		return make(ffparams, search, Parallelism.makeCpu(1));
	}
	
	public static MinimizingConfEnergyCalculator make(ForcefieldParams ffparams, SearchProblem search, Parallelism parallelism) {
		
		// make the forcefield interactions factory
		Factory<ForcefieldInteractions,Molecule> ffinteractions = (mol) -> FFInterGen.makeFullConf(search.confSpace, search.shellResidues, mol);
		
		// TODO: simplify this with a unified builder that uses the new Parallelism class
		// make the minimizer
		ConfMinimizer minimizer;
		switch (parallelism.type) {
			case Cpu:
				minimizer = new CpuConfMinimizer.Builder(ffparams, ffinteractions, search.confSpace)
					.setNumThreads(parallelism.numThreads)
					.build();
			break;
			case Gpu:
				minimizer = new GpuConfMinimizer.Builder(ffparams, ffinteractions, search.confSpace)
					.setGpuInfo(null, parallelism.numGpus, parallelism.numStreamsPerGpu)
					.build();
			break;
			default:
				throw new Error("unrecognized type: " + parallelism.type);
		}
		
		MinimizingConfEnergyCalculator ecalc = new MinimizingConfEnergyCalculator(minimizer);
		
		// add post pocessing steps
		if (search.useERef) {
			// the emat might not have been computed yet, so we can't get the ReferenceEnergies reference right now
			//ecalc.addConfPostProcessor(ConfPostProcessor.referenceEnergies(search.emat.geteRefMat()));
			// so look in the SearchProblem for the reference energies every time we post process a conf
			ecalc.addConfPostProcessor((econf) -> econf.offsetEnergy(-search.emat.geteRefMat().confERef(econf.getAssignments())));
		}
		if (search.addResEntropy) {
			ecalc.addConfPostProcessor(ConfPostProcessor.residueEntropy(search.confSpace));
		}
		
		return ecalc;
	}
	
	public static interface ConfPostProcessor {
		
		void postProcess(EnergiedConf conf);
		
		public static ConfPostProcessor referenceEnergies(ReferenceEnergies erefMat) {
			return (econf) -> econf.offsetEnergy(-erefMat.confERef(econf.getAssignments()));
		}
		
		public static ConfPostProcessor residueEntropy(ConfSpace confSpace) {
			return (econf) -> econf.offsetEnergy(confSpace.getConfResEntropy(econf.getAssignments()));
		}
	}
	
	private ConfMinimizer minimizer;
	private List<ConfPostProcessor> postProcessors;
	
	public MinimizingConfEnergyCalculator(ConfMinimizer minimizer) {
		this.minimizer = minimizer;
		this.postProcessors = new ArrayList<>();
	}
	
	public void addConfPostProcessor(ConfPostProcessor val) {
		this.postProcessors.add(val);
	}
	
	@Override
	public TaskExecutor getTasks() {
		return minimizer.getAsync().getTasks();
	}
	
	private EnergiedConf postProcessConf(EnergiedConf econf) {
		for (ConfPostProcessor postProcessor : postProcessors) {
			postProcessor.postProcess(econf);
		}
		return econf;
	}

	@Override
	public EnergiedConf calcEnergy(ScoredConf conf) {
		return postProcessConf(minimizer.getAsync().minimizeSync(conf));
	}
	
	@Override
	public void calcEnergyAsync(ScoredConf conf, Listener listener) {
		minimizer.getAsync().minimizeAsync(conf, new ConfMinimizer.Async.Listener() {
			@Override
			public void onMinimized(EnergiedConf econf) {
				listener.onFinished(postProcessConf(econf));
			}
		});
	}
	
	public void clean() {
		minimizer.cleanup();
	}
}
