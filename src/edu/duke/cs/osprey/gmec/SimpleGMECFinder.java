package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.externalMemory.EnergiedConfFIFOSerializer;
import edu.duke.cs.osprey.externalMemory.EnergiedConfPrioritySerializer;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.externalMemory.ScoredConfFIFOSerializer;
import edu.duke.cs.osprey.gmec.GMECFinder.ConfPruner;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.tools.*;

import java.io.*;
import java.util.*;
import java.util.function.Function;


/**
 * Searches a conformation space for the single conformation that minimizes the
 * desired energy function, ie the Global Minimum Energy Conformation, or GMEC.
 * 
 * The search is performed by ordering all conformations in the conformation space
 * using A* search, and evaluating the energy of each conformation in order.
 * 
 * The search terminates returns the lowest-energy conformation found so far
 * when it can be proven that all remaining conformations in the conformation
 * space must have higher energy.
 */
public class SimpleGMECFinder {

	public static final String ConfDBTableName = "GMEC";
	
	public static class Builder {
		
		/** A* implementation to sort conformations in the conformation space. */
		protected ConfSearch search;
		
		/** Calculates the energy for a conformation. */
		protected ConfEnergyCalculator confEcalc;
		protected ConfPruner pruner;
		protected ConfPrinter logPrinter;
		protected ConfPrinter consolePrinter;
		
		/** True to print all conformations found during A* search to the console */ 
		protected boolean printIntermediateConfsToConsole = false;
		
		/**
		 * True to use external memory (eg, disk, SSD, NAS) for when large data
		 * structures cannot fit in internal memory (eg, RAM).
		 * 
		 * Use {@link ExternalMemory#setInternalLimit} to set the amount of fixed internal memory
		 * and {@link ExternalMemory#setTempDir} to set the file path for external memory.
		 */
		protected boolean useExternalMemory = false;

		/**
		 * If a design experiences an unexpected abort, the conformation database can allow you to restore the
		 * design state and resume the calculation close to where it was aborted. Set a file to turn on the conf DB.
		 */
		protected File confDB = null;
		
		public Builder(ConfSearch search, ConfEnergyCalculator confEcalc) {
			this.search = search;
			this.confEcalc = confEcalc;
			this.pruner = null;
			this.logPrinter = new ConfPrinter.Nop();
			this.consolePrinter = new ConsoleConfPrinter();
		}
		
		/* TODO: integrate PartCR
		public Builder setConfPruner(ConfPruner val) {
			pruner = val;
			return this;
		}
		*/
		
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
		
		public Builder useExternalMemory() {
			ExternalMemory.checkInternalLimitSet();
			useExternalMemory = true;
			return this;
		}

		public Builder setConfDB(File val) {
			confDB = val;
			return this;
		}

		public SimpleGMECFinder build() {
			return new SimpleGMECFinder(
				search,
				confEcalc,
				pruner,
				logPrinter,
				consolePrinter,
				printIntermediateConfsToConsole,
				useExternalMemory,
				confDB
			);
		}
	}
	
	public ConfSearch search;
	public final ConfEnergyCalculator confEcalc;
	public final ConfPruner pruner;
	public final ConfPrinter logPrinter;
	public final ConfPrinter consolePrinter;
	public final boolean printIntermediateConfsToConsole;
	
	private final Queue.Factory.FIFO<ScoredConf> scoredFifoFactory;
	private final Queue.Factory.FIFO<EnergiedConf> energiedFifoFactory;
	private final Queue.Factory<EnergiedConf> energiedPriorityFactory;
	private final File confDBFile;

	protected SimpleGMECFinder(ConfSearch search, ConfEnergyCalculator confEcalc, ConfPruner pruner, ConfPrinter logPrinter, ConfPrinter consolePrinter, boolean printIntermediateConfsToConsole, boolean useExternalMemory, File confDBFile) {
		this.search = search;
		this.confEcalc = confEcalc;
		this.pruner = pruner;
		this.logPrinter = logPrinter;
		this.consolePrinter = consolePrinter;
		this.printIntermediateConfsToConsole = printIntermediateConfsToConsole;
		this.confDBFile = confDBFile;
		
		if (useExternalMemory) {
			RCs rcs = new RCs(confEcalc.confSpace);
			scoredFifoFactory = new Queue.ExternalFIFOFactory<>(new ScoredConfFIFOSerializer(rcs));
			energiedFifoFactory = new Queue.ExternalFIFOFactory<>(new EnergiedConfFIFOSerializer(rcs));
			energiedPriorityFactory = new Queue.ExternalPriorityFactory<>(new EnergiedConfPrioritySerializer(rcs));
		} else {
			scoredFifoFactory = new Queue.FIFOFactory<>();
			energiedFifoFactory = new Queue.FIFOFactory<>();
			energiedPriorityFactory = new Queue.PriorityFactory<>((a, b) -> Double.compare(a.getEnergy(), b.getEnergy()));
		}
	}

	private <T> T useDBIfNeeded(Function<ConfDB.ConfTable,T> block) {
		return ConfDB.useIfNeeded(confEcalc.confSpace, confDBFile, (confdb) -> {
			if (confdb == null) {
				return block.apply(null);
			} else {
				return block.apply(confdb.new ConfTable(ConfDBTableName));
			}
		});
	}

	public EnergiedConf find() {
		return find(0).poll();
	}
	
	public Queue.FIFO<EnergiedConf> find(double energyWindowSize) {

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

			// no confs in the search space
			System.out.println("No conformations found with finite energies");
			return Queue.FIFOFactory.of();
		}
		System.out.println("Found min score conformation in " + minScoreStopwatch.getTime(1));

		return useDBIfNeeded((confTable) -> {

			// evaluate the min score conf
			System.out.println("Computing energy of min score conf...");
			EnergiedConf eMinScoreConf = confEcalc.calcEnergy(minScoreConf, confTable);
			logPrinter.print(eMinScoreConf, confEcalc.confSpace);

			// peek ahead to the next conf
			ConfSearch.Splitter splitter = new ConfSearch.Splitter(search);
			ConfSearch.Splitter.Stream unpeekedConfs = splitter.makeStream();
			ConfSearch.Splitter.Stream peekedConfs = splitter.makeStream();
			ScoredConf peekedConf = peekedConfs.nextConf();
			peekedConfs.close();

			// do we need to check more confs?
			if (peekedConf == null) {

				// nope, there's no more confs, so we already have the GMEC
				System.out.println("Found GMEC! (it's actually the only conformation allowed by the conf space!)");
				consolePrinter.print(eMinScoreConf, confEcalc.confSpace);
				return Queue.FIFOFactory.of(eMinScoreConf);

			} else if (peekedConf.getScore() > eMinScoreConf.getEnergy() && energyWindowSize <= 0) {

				// nope, no confs have lower energy and we're not doing an energy window
				System.out.println("Found GMEC! (it's actually the min score conformation too!)");
				consolePrinter.print(eMinScoreConf, confEcalc.confSpace);
				return Queue.FIFOFactory.of(eMinScoreConf);

			} else {

				// yup, need to keep searching for the GMEC, or to enumerate an energy window
				// but drop the min score conf on the console before moving on
				consolePrinter.print(eMinScoreConf, confEcalc.confSpace);

				// estimate the top of our energy range
				/* NOTE:
					The "Energy Window" is the range of energies [GMEC energy,GMEC energy + window size].
					This range is merely an estimate of the true Energy Window based on the lowest energy we have so far.
					We'll refine this estimate as we evaluate more structures.
				*/
				final EnergyRange erange = new EnergyRange(eMinScoreConf.getEnergy(), energyWindowSize);

				// start the queue of energied confs
				Queue<EnergiedConf> econfs = energiedPriorityFactory.make();
				econfs.push(eMinScoreConf);

				checkMoreConfs(unpeekedConfs, erange, econfs, confTable);
				System.out.println(String.format("checked %d conformations", econfs.size()));

				// econfs are in a priority queue, so the first one is the GMEC
				EnergiedConf gmec = econfs.peek();
				System.out.println("\nFound GMEC!");
				consolePrinter.print(gmec, confEcalc.confSpace);

				// return just the confs in the energy window
				Queue.FIFO<EnergiedConf> econfsInRange = econfs.filterTo(
					energiedFifoFactory.make(),
					(conf) -> erange.contains(conf.getEnergy())
				);

				System.out.println(String.format("Found %d total conformations in energy window", econfsInRange.size()));

				return econfsInRange;
			}
		});
	}
	
	private void checkMoreConfs(ConfSearch search, EnergyRange erange, Queue<EnergiedConf> econfs, ConfDB.ConfTable confTable) {
		
		setErangeProgress(search, erange);
		
		// enumerate all confs in order of the scores, up to the estimate of the top of the energy window
		System.out.println("Enumerating other low-scoring conformations...");
		Queue.FIFO<ScoredConf> otherLowEnergyConfs = scoredFifoFactory.make();
		Stopwatch timingStopwatch = new Stopwatch().start();
		Stopwatch speculativeMinimizationStopwatch = new Stopwatch().start();
		while (true) {
			
			// get the next conf, or stop searching if none left
			ScoredConf conf = search.nextConf();
			if (conf == null) {
				break;
			}
			
			// ignore the conf if it's out of range
			if (conf.getScore() > erange.getMax()) {
				break;
			}
			
			// save the conf for later minimization
			otherLowEnergyConfs.push(conf);
			
			// if we're exactly at the limit, stop after saving the conf
			if (conf.getScore() == erange.getMax()) {
				break;
			}
			
			// if we've been enumerating confs for a while, try a minimization to see if we get a smaller window
			if (speculativeMinimizationStopwatch.getTimeS() >= 10) {
				speculativeMinimizationStopwatch.stop();
				
				// save the conf and the energy for later
				EnergiedConf econf = confEcalc.calcEnergy(otherLowEnergyConfs.poll(), confTable);
				handleEnergiedConf(econf, econfs, erange);
				
				boolean changed = erange.updateMin(econf.getEnergy());
				if (changed) {
					System.out.println(String.format("Lower conformation energy updated energy window! remaining: %14.8f", erange.getMax() - conf.getScore()));
					setErangeProgress(search, erange);
				}
				
				speculativeMinimizationStopwatch.start();
			}
		}
		
		System.out.println(String.format("\tFound %d more in %s", otherLowEnergyConfs.size(), timingStopwatch.getTime(1)));
		
		if (!otherLowEnergyConfs.isEmpty()) {
			
			/* TODO: integrate PartCR
			if (pruner != null) {
				pruner.prune(otherLowEnergyConfs, ecalc);
			}
			*/
			
			minimizeLowEnergyConfs(otherLowEnergyConfs, erange, econfs, confTable);
		}
	}

	private void minimizeLowEnergyConfs(Queue.FIFO<ScoredConf> lowEnergyConfs, EnergyRange erange, Queue<EnergiedConf> econfs, ConfDB.ConfTable confTable) {

		// calculate energy for each conf
		// this will probably take a while, so track progress
		Progress progress = new Progress(lowEnergyConfs.size());

		// what to do when we get a conf energy?
		TaskListener<EnergiedConf> ecalcListener = (econf) -> {

			// NOTE: this is called on a listener thread, which is separate from the main thread

			handleEnergiedConf(econf, econfs, erange);

			// refine the estimate of the top of the energy window
			boolean changed = erange.updateMin(econf.getEnergy());
			if (changed) {

				System.out.println(String.format("\nNew lowest energy: %.6f", erange.getMin()));

				// remove confs from the queue whose energies are above the new lower window
				// but don't race the main thread which is poll()ing the queue
				synchronized (lowEnergyConfs) {

					long numLowEnergyConfsBefore = lowEnergyConfs.size();

					// prune conformations with the new window
					// (in-place, since this is a FIFO queue)
					lowEnergyConfs.filter((conf) -> erange.containsOrBelow(conf.getScore()));

					if (lowEnergyConfs.size() < numLowEnergyConfsBefore) {
						System.out.println(String.format("\tReduced to %d low-energy conformations", lowEnergyConfs.size()));
					}
					progress.setTotalWork(lowEnergyConfs.size());
				}
			}

			progress.incrementProgress();
		};

		// calc the conf energy asynchronously
		System.out.println(String.format("\nComputing energies for %d conformations...", lowEnergyConfs.size()));
		while (true) {

			// get the next conf to calc the energy for
			// but don't race the listener thread, which can sometimes filter the queue
			ScoredConf conf;
			synchronized (lowEnergyConfs) {

				if (lowEnergyConfs.isEmpty()) {
					break;
				}

				conf = lowEnergyConfs.poll();
			}

			// send the conf to the energy calculator
			confEcalc.calcEnergyAsync(conf, confTable, ecalcListener);
		}

		confEcalc.tasks.waitForFinish();
	}

	private void setErangeProgress(ConfSearch confSearch, EnergyRange erange) {
		
		// HACKHACK: set progress goal
		if (confSearch instanceof ConfAStarTree) {
			ConfAStarTree tree = (ConfAStarTree)confSearch;
			if (tree.getProgress() != null) {
				tree.getProgress().setGoalScore(erange.getMax());
			}
		} else if (confSearch instanceof ConfSearch.Splitter.Stream) {
			setErangeProgress(((ConfSearch.Splitter.Stream)confSearch).getSource(), erange);
		}
	}
	
	private void handleEnergiedConf(EnergiedConf econf, Queue<EnergiedConf> econfs, EnergyRange erange) {
		
		// make sure the score was actually a lower bound
		if (econf.getScore() > econf.getEnergy() + 0.1) {
			System.out.println(String.format("WARNING: Conformation score (%f) is not a lower bound on the energy (%f)."
				+ "\n\tThis is evidence that OSPREY is unable to guarantee finding exactly the GMEC for this design,"
				+ " but we'll probably get it anyway, or at least get really close."
				+ "\n\tAssignments: %s",
				econf.getScore(),
				econf.getEnergy(),
				Arrays.toString(econf.getAssignments())
			));
		}

		// save the conf and the energy for later
		econfs.push(econf);

		// immediately output the conf, in case the run aborts and we want to resume later
		logPrinter.print(econf, confEcalc.confSpace);

		// log the conf to console too if desired
		if (printIntermediateConfsToConsole) {
			System.out.println();
			consolePrinter.print(econf, confEcalc.confSpace, erange);
		}
	}
}
