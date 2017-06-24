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
