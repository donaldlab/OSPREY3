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

import static edu.duke.cs.osprey.tools.Log.formatBig;


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

		/** Print status and progress info to the console */
		protected boolean printToConsole = true;
		
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

		public Builder setPrintToConsole(boolean val) {
			printToConsole = val;
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
				printToConsole,
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
	public final boolean printToConsole;
	
	private final Queue.Factory.FIFO<ScoredConf> scoredFifoFactory;
	private final Queue.Factory.FIFO<EnergiedConf> energiedFifoFactory;
	private final Queue.Factory<EnergiedConf> energiedPriorityFactory;
	private final File confDBFile;

	protected SimpleGMECFinder(ConfSearch search, ConfEnergyCalculator confEcalc, ConfPruner pruner, ConfPrinter logPrinter, ConfPrinter consolePrinter, boolean printIntermediateConfsToConsole, boolean printToConsole, boolean useExternalMemory, File confDBFile) {
		this.search = search;
		this.confEcalc = confEcalc;
		this.pruner = pruner;
		this.logPrinter = logPrinter;
		this.consolePrinter = consolePrinter;
		this.printIntermediateConfsToConsole = printIntermediateConfsToConsole;
		this.printToConsole = printToConsole;
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

	private void log(String msg, Object ... args) {
		if (printToConsole) {
			edu.duke.cs.osprey.tools.Log.log(msg, args);
		}
	}

	public EnergiedConf find() {
		return find(0).poll();
	}
	
	public Queue.FIFO<EnergiedConf> find(double energyWindowSize) {

		// start searching for the min score conf
		log("Searching for min score conformation...");
		Stopwatch minScoreStopwatch = new Stopwatch().start();
		try {
			log("\t(among %s possibilities)", formatBig(search.getNumConformations()));
		} catch (UnsupportedOperationException ex) {
			// conf search doesn't support it, no big deal
		}
		ScoredConf minScoreConf = search.nextConf();
		if (minScoreConf == null) {

			// no unpruned confs in the search space
			log("No conformations found with finite energies (possibly all confs have been pruned)");
			return Queue.FIFOFactory.of();
		}
		log("Found min score conformation in %s", minScoreStopwatch.getTime(1));

		// open the ConfDB if needed
		try (ConfDB confdb = ConfDB.makeIfNeeded(confEcalc.confSpace, confDBFile)) {
			ConfDB.ConfTable confTable = null;
			if (confdb != null) {
				confTable = confdb.new ConfTable(ConfDBTableName);
			}

			// evaluate the min score conf
			log("Computing energy of min score conf...");
			EnergiedConf eMinScoreConf = confEcalc.calcEnergy(minScoreConf, confTable);
			logPrinter.print(eMinScoreConf, confEcalc.confSpace);

			// peek ahead to the next conf
			ConfSearch.MultiSplitter splitter = new ConfSearch.MultiSplitter(search);
			ConfSearch.MultiSplitter.Stream unpeekedConfs = splitter.makeStream();
			ConfSearch.MultiSplitter.Stream peekedConfs = splitter.makeStream();
			ScoredConf peekedConf = peekedConfs.nextConf();
			peekedConfs.close();

			// do we need to check more confs?
			if (peekedConf == null) {

				// nope, there's no more confs, so we already have the GMEC
				log("Found GMEC! (it's actually the only conformation allowed by the conf space!)");
				if (printToConsole) {
					consolePrinter.print(eMinScoreConf, confEcalc.confSpace);
				}
				return Queue.FIFOFactory.of(eMinScoreConf);

			} else if (peekedConf.getScore() > eMinScoreConf.getEnergy() && energyWindowSize <= 0) {

				// nope, no confs have lower energy and we're not doing an energy window
				log("Found GMEC! (it's actually the min score conformation too!)");
				if (printToConsole) {
					consolePrinter.print(eMinScoreConf, confEcalc.confSpace);
				}
				return Queue.FIFOFactory.of(eMinScoreConf);

			} else {

				// yup, need to keep searching for the GMEC, or to enumerate an energy window
				// but drop the min score conf on the console before moving on
				if (printToConsole) {
					consolePrinter.print(eMinScoreConf, confEcalc.confSpace);
				}

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
				log("checked %d conformations", econfs.size());

				// econfs are in a priority queue, so the first one is the GMEC
				EnergiedConf gmec = econfs.peek();
				log("\nFound GMEC!");
				if (printToConsole) {
					consolePrinter.print(gmec, confEcalc.confSpace);
				}

				// return just the confs in the energy window
				Queue.FIFO<EnergiedConf> econfsInRange = econfs.filterTo(
					energiedFifoFactory.make(),
					(conf) -> erange.contains(conf.getEnergy())
				);

				log("Found %d total conformations in energy window", econfsInRange.size());

				return econfsInRange;
			}
		}
	}
	
	private void checkMoreConfs(ConfSearch search, EnergyRange erange, Queue<EnergiedConf> econfs, ConfDB.ConfTable confTable) {
		
		setErangeProgress(search, erange);
		
		// enumerate all confs in order of the scores, up to the estimate of the top of the energy window
		log("Enumerating other low-scoring conformations...");
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
					log("Lower conformation energy updated energy window! remaining: %14.8f", erange.getMax() - conf.getScore());
					setErangeProgress(search, erange);
				}
				
				speculativeMinimizationStopwatch.start();
			}
		}
		
		log("\tFound %d more in %s", otherLowEnergyConfs.size(), timingStopwatch.getTime(1));
		
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
		Progress progress;
		if (printToConsole) {
			progress = new Progress(lowEnergyConfs.size());
		} else {
			progress = null;
		}

		// what to do when we get a conf energy?
		TaskListener<EnergiedConf> ecalcListener = (econf) -> {

			// NOTE: this is called on a listener thread, which is separate from the main thread

			handleEnergiedConf(econf, econfs, erange);

			// refine the estimate of the top of the energy window
			boolean changed = erange.updateMin(econf.getEnergy());
			if (changed) {

				log("\nNew lowest energy: %.6f", erange.getMin());

				// remove confs from the queue whose energies are above the new lower window
				// but don't race the main thread which is poll()ing the queue
				synchronized (lowEnergyConfs) {

					long numLowEnergyConfsBefore = lowEnergyConfs.size();

					// prune conformations with the new window
					// (in-place, since this is a FIFO queue)
					lowEnergyConfs.filter((conf) -> erange.containsOrBelow(conf.getScore()));

					if (lowEnergyConfs.size() < numLowEnergyConfsBefore) {
						log("\tReduced to %d low-energy conformations", lowEnergyConfs.size());
					}
					if (progress != null) {
						progress.setTotalWork(lowEnergyConfs.size());
					}
				}
			}

			if (progress != null) {
				progress.incrementProgress();
			}
		};

		// calc the conf energy asynchronously
		log("\nComputing energies for %d conformations...", lowEnergyConfs.size());
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
		} else if (confSearch instanceof ConfSearch.MultiSplitter.Stream) {
			setErangeProgress(((ConfSearch.MultiSplitter.Stream)confSearch).getSource(), erange);
		}
	}
	
	private void handleEnergiedConf(EnergiedConf econf, Queue<EnergiedConf> econfs, EnergyRange erange) {
		
		// make sure the score was actually a lower bound
		if (econf.getScore() > econf.getEnergy() + 0.1) {
			log("WARNING: Conformation score (%f) is not a lower bound on the energy (%f)."
				+ "\n\tThis is evidence that OSPREY is unable to guarantee finding exactly the GMEC for this design,"
				+ " but we'll probably get it anyway, or at least get really close."
				+ "\n\tAssignments: %s",
				econf.getScore(),
				econf.getEnergy(),
				Arrays.toString(econf.getAssignments())
			);
		}

		// save the conf and the energy for later
		econfs.push(econf);

		// immediately output the conf, in case the run aborts and we want to resume later
		logPrinter.print(econf, confEcalc.confSpace);

		// log the conf to console too if desired
		if (printToConsole && printIntermediateConfsToConsole) {
			log("");
			consolePrinter.print(econf, confEcalc.confSpace, erange);
		}
	}
}
