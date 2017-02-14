package edu.duke.cs.osprey.gmec;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.gmec.GMECFinder.ConfPruner;
import edu.duke.cs.osprey.minimization.SimpleConfMinimizer;
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
		
		public Builder setPrintIntermediateConfsToConsole(boolean val) {
			printIntermediateConfsToConsole = val;
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
	
	private SimpleGMECFinder(SimpleConfSpace space, ConfSearch search, ConfEnergyCalculator.Async ecalc, ConfPruner pruner, ConfPrinter logPrinter, ConfPrinter consolePrinter, boolean printIntermediateConfsToConsole) {
		this.space = space;
		this.search = search;
		this.ecalc = ecalc;
		this.pruner = pruner;
		this.logPrinter = logPrinter;
		this.consolePrinter = consolePrinter;
		this.printIntermediateConfsToConsole = printIntermediateConfsToConsole;
	}
	
	public EnergiedConf find() {
		
		List<EnergiedConf> confs = find(0);
		
		if (!confs.isEmpty()) {
			return confs.get(0);
		}
		
		return null;
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
		System.out.println("Computing energy of min score conf...");
		EnergiedConf eMinScoreConf = ecalc.calcEnergy(minScoreConf);
		logPrinter.print(eMinScoreConf, space);
		
		// peek ahead to the next conf
		ConfSearch.Splitter splitter = new ConfSearch.Splitter(search);
		ConfSearch unpeekedConfs = splitter.makeStream();
		ScoredConf peekedConf = splitter.makeStream().nextConf();
		
		// do we need to check more confs?
		if (peekedConf == null) {
			
			// nope, there's no more confs, so we already have the GMEC
			System.out.println("Found GMEC! (it's actually the only conformation allowed by the conf space!)");
			consolePrinter.print(eMinScoreConf, space);
			return Arrays.asList(eMinScoreConf);
			
		} else if (peekedConf.getScore() > eMinScoreConf.getEnergy() && energyWindowSize <= 0) {
			
			// nope, no confs have lower energy and we're not doing an energy window
			System.out.println("Found GMEC!");
			consolePrinter.print(eMinScoreConf, space);
			return Arrays.asList(eMinScoreConf);
		
		} else {
			
			// yup, need to keep searching for the GMEC, or to enumerate an energy window
			// but drop the min score conf on the console before moving on
			consolePrinter.print(eMinScoreConf, space);
			
			// estimate the top of our energy range
			// this is an upper bound for now, we'll refine it as we evaluate more structures
			final EnergyRange erange = new EnergyRange(eMinScoreConf.getEnergy(), energyWindowSize);
		
			// start the list of energied confs
			List<EnergiedConf> econfs = new ArrayList<>();
			econfs.add(eMinScoreConf);
		
			checkMoreConfs(unpeekedConfs, erange, econfs);
			
			// econfs is sorted now, so the first one is the GMEC
			EnergiedConf gmec = econfs.get(0);
			System.out.println("\nFound GMEC!");
			consolePrinter.print(gmec, space);
			
			System.out.println(String.format("Found %d total conformations in energy window", econfs.size()));
			
			return econfs;
		}
	}
	
	private void checkMoreConfs(ConfSearch search, EnergyRange erange, List<EnergiedConf> econfs) {
		
		// enumerate all confs in order of the scores, up to the estimate of the top of the energy window
		System.out.println("Enumerating other low-scoring conformations...");
		List<ScoredConf> otherLowEnergyConfs = new ArrayList<>(search.nextConfs(erange.getMax()));
		
		// ConfSearch implementations can give one extra conf that's technically above the bound we gave to nextConfs()
		// that's because all ConfSearch implementations aren't smart enough (yet) to peek ahead
		// so prune the conf now because we definitely don't care about it, or any other higher confs
		pruneConfsOutsideRange(otherLowEnergyConfs, erange);
		
		System.out.println(String.format("\tFound %d more", otherLowEnergyConfs.size()));
		
		if (!otherLowEnergyConfs.isEmpty()) {
			
			// prune the confs list if needed (eg, PartCR)
			if (pruner != null) {
				pruner.prune(otherLowEnergyConfs, ecalc);
			}

			minimizeLowEnergyConfs(otherLowEnergyConfs, erange, econfs);
		}

		// sort all the confs by energy
		Collections.sort(econfs, new Comparator<EnergiedConf>() {
			@Override
			public int compare(EnergiedConf a, EnergiedConf b) {
				return Double.compare(a.getEnergy(), b.getEnergy());
			}
		});
		
		// prune all confs outside the energy window
		pruneConfsOutsideRange(econfs, erange);
	}

	private void minimizeLowEnergyConfs(List<ScoredConf> lowEnergyConfs, EnergyRange erange, List<EnergiedConf> econfs) {
		
		// calculate energy for each conf
		// this will probably take a while, so track progress
		Progress progress = new Progress(lowEnergyConfs.size());

		// what to do when we get a conf energy?
		ConfEnergyCalculator.Async.Listener ecalcListener = (econf) -> {
			
			// make sure the score was actually a lower bound
			if (econf.getScore() > econf.getEnergy() + 0.1) {
				throw new Error(String.format("Conformation score (%f) is not a lower bound on the energy (%f)! This is a serious bug.",
					econf.getScore(),
					econf.getEnergy()
				));
			}

			// save the conf and the energy for later
			econfs.add(econf);

			// immediately output the conf, in case the run aborts and we want to resume later
			System.out.println();
			logPrinter.print(econf, space);

			// log the conf to console too if desired
			if (printIntermediateConfsToConsole) {
				consolePrinter.print(econf, space, erange);
			}

			progress.incrementProgress();

			// refine the estimate of the top of the energy window
			boolean changed = erange.updateMin(econf.getEnergy());
			if (changed) {

				// prune conformations with the new window
				pruneConfsOutsideRange(lowEnergyConfs, erange);

				// update progress
				System.out.println(String.format("\nNew lowest energy: %.6f", erange.getMin()));
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
	
	private void pruneConfsOutsideRange(List<? extends ScoredConf> confs, EnergyRange erange) {
		for (int i=confs.size()-1; i>=0; i--) {
			if (confs.get(i).getScore() > erange.getMax()) {
				confs.remove(i);
			}
		}
	}
}
