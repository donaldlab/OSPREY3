package edu.duke.cs.osprey.energy;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.Progress;

public class MinimizingConfEnergyCalculator implements ConfEnergyCalculator.Async {
	
	public static class Builder {
		
		private FragmentEnergyCalculator.Async ecalc;
		private SimpleReferenceEnergies eref = null;
		
		public Builder(FragmentEnergyCalculator.Async ecalc) {
			this.ecalc = ecalc;
		}
		
		public Builder setReferenceEnegries(SimpleReferenceEnergies val) {
			this.eref = val;
			return this;
		}
		
		public MinimizingConfEnergyCalculator build() {
			return new MinimizingConfEnergyCalculator(ecalc, eref);
		}
	}
	
	// TODO: entropy energies
	
	public final FragmentEnergyCalculator.Async ecalc;
	public final SimpleReferenceEnergies eref;
	
	private MinimizingConfEnergyCalculator(FragmentEnergyCalculator.Async ecalc, SimpleReferenceEnergies eref) {
		this.ecalc = ecalc;
		this.eref = eref;
	}
	
	@Override
	public EnergiedConf calcEnergy(ScoredConf conf) {
		
		RCTuple frag = new RCTuple(conf.getAssignments());
		
		ResInterGen gen = ResInterGen.of(ecalc.getConfSpace())
			.addIntras(frag)
			.addInters(frag)
			.addShell(frag);
		
		// add reference energies if needed
		if (eref != null) {
			gen.add(eref.getFragmentEnergy(ecalc.getConfSpace(), frag));
		}
		
		double energy = ecalc.calcEnergy(frag, gen.make());
		
		return new EnergiedConf(conf, energy);
	}

	@Override
	public void calcEnergyAsync(ScoredConf conf, ConfEnergyCalculator.Async.Listener listener) {
		ecalc.getTasks().submit(() -> calcEnergy(conf), listener);
	}
	
	@Override
	public TaskExecutor getTasks() {
		return ecalc.getTasks();
	}

	@Override
	public void cleanup() {
		ecalc.cleanup();
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
		getTasks().waitForFinish();
		
		return econfs;
	}
}
