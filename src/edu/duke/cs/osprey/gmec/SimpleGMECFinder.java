package edu.duke.cs.osprey.gmec;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.gmec.GMECFinder.ConfPruner;
import edu.duke.cs.osprey.minimization.SimpleConfMinimizer;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Progress;
import edu.duke.cs.osprey.tools.Stopwatch;

public class SimpleGMECFinder {
	
	public static class Builder {
		
		private SimpleConfSpace space;
		private ConfSearch search;
		private ConfEnergyCalculator.Async ecalc;
		private ConfPruner pruner;
		private ConfPrinter logPrinter;
		private ConfPrinter consolePrinter;
		private boolean printIntermediateConfsToConsole;
		
		public Builder(SimpleConfSpace space, ConfSearch search) {
			this.space = space;
			this.search = search;
			this.ecalc = SimpleConfMinimizer.builder(space).build();
			this.pruner = null;
			this.logPrinter = new ConfPrinter.Nop();
			this.consolePrinter = new ConsoleConfPrinter();
			this.printIntermediateConfsToConsole = false;
		}
		
		public Builder setEnergyCalculator(ConfEnergyCalculator val) {
			return setEnergyCalculator(new ConfEnergyCalculator.Async.Adapter(val));
		}
		
		public Builder setEnergyCalculator(ConfEnergyCalculator.Async val) {
			ecalc = val;
			return this;
		}
		
		public Builder setConfPruner(ConfPruner val) {
			pruner = val;
			return this;
		}
		
		public Builder setLogPrinter(ConfPrinter val) {
			logPrinter = val;
			return this;
		}
		
		public Builder setConsolePrinter(ConfPrinter val) {
			consolePrinter = val;
			return this;
		}
		
		public SimpleGMECFinder build() {
			return new SimpleGMECFinder(
				space,
				search,
				ecalc,
				pruner,
				logPrinter,
				consolePrinter,
				printIntermediateConfsToConsole
			);
		}
	}
	
	public static Builder builder(SimpleConfSpace confSpace, EnergyMatrix emat) {
		return new Builder(confSpace, ConfAStarTree.builder(emat, confSpace).build());
	}

	public static Builder builder(SimpleConfSpace confSpace, ConfSearch confSearch) {
		return new Builder(confSpace, confSearch);
	}
	
	public final SimpleConfSpace space;
	public final ConfSearch search;
	public final ConfEnergyCalculator.Async ecalc;
	public final ConfPruner pruner;
	public final ConfPrinter logPrinter;
	public final ConfPrinter consolePrinter;
	public final boolean printIntermediateConfsToConsole;
	
	public Molecule molMinGMEC;
	
	private SimpleGMECFinder(SimpleConfSpace space, ConfSearch search, ConfEnergyCalculator.Async ecalc, ConfPruner pruner, ConfPrinter logPrinter, ConfPrinter consolePrinter, boolean printIntermediateConfsToConsole) {
		
		this.space = space;
		this.search = search;
		this.ecalc = ecalc;
		this.pruner = pruner;
		this.logPrinter = logPrinter;
		this.consolePrinter = consolePrinter;
		this.printIntermediateConfsToConsole = printIntermediateConfsToConsole;
		
		molMinGMEC = null;
	}
	
	public EnergiedConf find() {
		
		List<EnergiedConf> confs = find(0);
		if (confs.isEmpty()) {
			return null;
		}
		
		return confs.get(0);
	}
	
	public List<EnergiedConf> find(double energyWindowSize) {
		
		// start searching for the min score conf
		System.out.println("Searching for min score conformation...");
		Stopwatch minScoreStopwatch = new Stopwatch().start();
		try {
			System.out.println("\t(among " + search.getNumConformations().floatValue() + " possibilities)");
		} catch (UnsupportedOperationException ex) {
			// conf search doesn't support it, no big deal
		}
		ScoredConf minScoreConf = search.nextConf();
		if (minScoreConf == null) {
			
			// no confs in the search space, can't recover, just bail
			return new ArrayList<>();
		}
		System.out.println("Found min score conformation in " + minScoreStopwatch.getTime(1));
		
		// evaluate the min score conf
		System.out.println("Computing energy...");
		EnergiedConf eMinScoreConf = ecalc.calcEnergy(minScoreConf);
		logPrinter.print(space, eMinScoreConf);
		consolePrinter.print(space, eMinScoreConf);
		
		List<EnergiedConf> econfs = new ArrayList<>();
		econfs.add(eMinScoreConf);
		
		// estimate the top of our energy window
		// this is an upper bound for now, we'll refine it as we evaluate more structures
		final EnergyWindow window = new EnergyWindow(eMinScoreConf.getEnergy(), energyWindowSize);
		
		// enumerate all confs in order of the scores, up to the estimate of the top of the energy window
		System.out.println("Enumerating other low-scoring conformations...");
		List<ScoredConf> lowEnergyConfs = new ArrayList<>();
		lowEnergyConfs.add(minScoreConf);
		lowEnergyConfs.addAll(search.nextConfs(window.getMax()));
		System.out.println(String.format("\tFound %d more", lowEnergyConfs.size() - 1));
		
		if (!lowEnergyConfs.isEmpty()) {

			// prune the confs list
			if (pruner != null) {
				pruner.prune(lowEnergyConfs, ecalc);
			}

			// calculate energy for each conf
			// this will probably take a while, so track progress
			Progress progress = new Progress(lowEnergyConfs.size());

			// what to do when we get a conf energy?
			ConfEnergyCalculator.Async.Listener ecalcListener = (econf) -> {

				// save the conf and the energy for later
				econfs.add(econf);

				// immediately output the conf, in case the run aborts and we want to resume later
				logPrinter.print(space, econf);

				// log the conf to console if desired
				if (printIntermediateConfsToConsole) {
					System.out.println("\nENUMERATING CONFORMATION");
					consolePrinter.print(space, econf, window);
				}

				progress.incrementProgress();

				// refine the estimate of the top of the energy window
				boolean changed = window.update(econf.getEnergy());
				if (changed) {

					// prune conformations with the new window
					for (int i=lowEnergyConfs.size()-1; i>=0; i--) {
						if (lowEnergyConfs.get(i).getScore() > window.getMax()) {
							lowEnergyConfs.remove(i);
						} else {
							break;
						}
					}

					// update progress
					System.out.println(String.format("\nNew lowest energy: %.6f", window.getMin()));
					System.out.println(String.format("\tReduced to %d low-energy conformations", lowEnergyConfs.size()));
					progress.setTotalWork(lowEnergyConfs.size());
				}
			};

			// calc the conf energy asynchronously
			System.out.println(String.format("\nComputing energies for %d conformations...", lowEnergyConfs.size()));
			for (int i=0; i<lowEnergyConfs.size(); i++) {
				ecalc.calcEnergyAsync(lowEnergyConfs.get(i), ecalcListener);
			}
			ecalc.waitForFinish();
		}
		
		// sort all the confs by energy
		Collections.sort(econfs, new Comparator<EnergiedConf>() {
			@Override
			public int compare(EnergiedConf a, EnergiedConf b) {
				return Double.compare(a.getEnergy(), b.getEnergy());
			}
		});
		
		// minEnergyConf is the minGMEC!! =)
		EnergiedConf minGMEC = econfs.get(0);
		System.out.println("\nFound minGMEC!");
		consolePrinter.print(space, minGMEC);
		
		assert (minGMEC.getEnergy() == window.getMin());
		
		// prune all confs outside the energy window and return them
		Iterator<EnergiedConf> iter = econfs.iterator();
		while (iter.hasNext()) {
			EnergiedConf econf = iter.next();
			if (!window.contains(econf.getEnergy())) {
				iter.remove();
			}
		}
		
		if (econfs.size() > 1) {
			System.out.println(String.format("Also found %d more conformations in energy window", econfs.size() - 1));
		}
		
		return econfs;
	}
}
