package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.ForcefieldInteractionsGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.minimization.ConfMinimizer;
import edu.duke.cs.osprey.minimization.CpuConfMinimizer;
import edu.duke.cs.osprey.minimization.GpuConfMinimizer;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public class MinimizingEnergyCalculator implements ConfEnergyCalculator.Async {
	
	// TODO: this should eventually go into a CFP-only area
	// it can be moved when we start refactoring config stuff to prepare for Python-land
	public static MinimizingEnergyCalculator makeFromConfig(SearchProblem search, ConfigFileParser cfp, boolean areConfsStreamed) {
		int numThreads = cfp.getParams().getInt("MinimizationThreads");
		int numGpus = cfp.getParams().getInt("MinimizationGpus");
		int streamsPerGpu = cfp.getParams().getInt("MinimizationStreamsPerGpu");
		return make(EnvironmentVars.curEFcnGenerator.ffParams, search, numGpus, streamsPerGpu, numThreads, areConfsStreamed);
	}
	
	public static MinimizingEnergyCalculator make(ForcefieldParams ffparams, SearchProblem search, int numGpus, int streamsPerGpu, int numThreads, boolean areConfsStreaming) {
		
		// make the forcefield interactions factory
		ForcefieldInteractionsGenerator intergen = new ForcefieldInteractionsGenerator();
		Factory<ForcefieldInteractions,Molecule> ffinteractions = (mol) -> intergen.makeFullConf(search.confSpace, search.shellResidues, mol);
		
		ConfMinimizer minimizer;
		if (numGpus > 0) {
			minimizer = new GpuConfMinimizer.Builder(ffparams, ffinteractions, search.confSpace)
				.setGpuInfo(null, numGpus, streamsPerGpu)
				.setAreConfsStreaming(areConfsStreaming)
				.build();
		} else {
			minimizer = new CpuConfMinimizer.Builder(ffparams, ffinteractions, search.confSpace)
				.setAreConfsStreaming(areConfsStreaming)
				.setNumThreads(numThreads)
				.build();
		}
		
		return new MinimizingEnergyCalculator(search, minimizer);
	}
		
	private SearchProblem search;
	private ConfMinimizer minimizer;
	
	public MinimizingEnergyCalculator(SearchProblem search, ConfMinimizer minimizer) {
		this.search = search;
		this.minimizer = minimizer;
	}
	
	@Override
	public int getParallelism() {
		return minimizer.getAsync().getParallelism();
	}
	
	private EnergiedConf postProcessConf(EnergiedConf econf) {
		
		// add post-minimization energy modifications
		if (search.useERef) {
			econf.offsetEnergy(-search.emat.geteRefMat().confERef(econf.getAssignments()));
		}
		if (search.addResEntropy) {
			econf.offsetEnergy(search.confSpace.getConfResEntropy(econf.getAssignments()));
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
				listener.onEnergy(postProcessConf(econf));
			}
		});
	}
	
	@Override
	public void waitForSpace() {
		minimizer.getAsync().waitForSpace();
	}
	
	@Override
	public void waitForFinish() {
		minimizer.getAsync().waitForFinish();
	}

	@Override
	public void cleanup() {
		minimizer.cleanup();
	}
}
