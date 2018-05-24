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

package edu.duke.cs.osprey.lute;


import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.externalMemory.*;
import edu.duke.cs.osprey.gmec.ConfPrinter;
import edu.duke.cs.osprey.gmec.ConsoleConfPrinter;
import edu.duke.cs.osprey.gmec.EnergyRange;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;


/**
 * Just like SimpleGMECFinder, except it doesn't use lower bounds on conformation
 * energies at all. Instead, this GMEC finder assumes LUTE energies are the minimized
 * energies, and enumerates the confs in order of LUTE energies using A* search.
 */
public class LUTEGMECFinder {

	public static class Builder {

		private PruningMatrix pmat;

		/** Calculates the energy for a conformation using the information learned by LUTE. */
		private LUTEConfEnergyCalculator confEcalc;

		private ConfPrinter logPrinter;
		private ConfPrinter consolePrinter;

		/** True to print all conformations found during A* search to the console */
		private boolean printIntermediateConfsToConsole = false;

		public Builder(PruningMatrix pmat, LUTEConfEnergyCalculator confEcalc) {
			this.pmat = pmat;
			this.confEcalc = confEcalc;
			this.logPrinter = new ConfPrinter.Nop();
			this.consolePrinter = new ConsoleConfPrinter();
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

		public LUTEGMECFinder build() {
			return new LUTEGMECFinder(
				pmat,
				confEcalc,
				logPrinter,
				consolePrinter,
				printIntermediateConfsToConsole
			);
		}
	}

	public final PruningMatrix pmat;
	public final LUTEConfEnergyCalculator confEcalc;
	public final ConfPrinter logPrinter;
	public final ConfPrinter consolePrinter;
	public final boolean printIntermediateConfsToConsole;

	private LUTEGMECFinder(PruningMatrix pmat, LUTEConfEnergyCalculator confEcalc, ConfPrinter logPrinter, ConfPrinter consolePrinter, boolean printIntermediateConfsToConsole) {
		this.pmat = pmat;
		this.confEcalc = confEcalc;
		this.logPrinter = logPrinter;
		this.consolePrinter = consolePrinter;
		this.printIntermediateConfsToConsole = printIntermediateConfsToConsole;
	}

	public ConfSearch.ScoredConf find() {
		return find(0).poll();
	}

	public Queue.FIFO<ConfSearch.ScoredConf> find(double energyWindowSize) {

		Queue.FIFO<ConfSearch.ScoredConf> confs = Queue.FIFOFactory.of();

		ConfAStarTree search = new ConfAStarTree.Builder(null, pmat)
			.setLUTE(confEcalc)
			.build();

		// start searching for the min score conf
		System.out.println("Searching for GMEC...");
		Stopwatch minScoreStopwatch = new Stopwatch().start();
		System.out.println("\t(among " + search.getNumConformations().floatValue() + " possibilities)");
		ConfSearch.ScoredConf gmec = search.nextConf();
		if (gmec == null) {

			// no confs in the search space
			System.out.println("No conformations found with finite energies");
			return confs;
		}
		System.out.println("Found GMEC in " + minScoreStopwatch.getTime(1));
		consolePrinter.print(gmec, confEcalc.confSpace);
		logPrinter.print(gmec, confEcalc.confSpace);
		confs.push(gmec);

		// do we need to enumate the energy window?
		if (energyWindowSize > 0) {

			// yup, make the energy window
			final EnergyRange ewin = new EnergyRange(gmec.getScore(), energyWindowSize);

			while (true) {

				ConfSearch.ScoredConf conf = search.nextConf();
				if (conf == null) {

					// no confs left in the search space
					System.out.println("Conformation space exhausted, can't finished enumerating energy window.");
					break;
				}

				if (ewin.contains(conf.getScore())) {

					// record the conf
					confs.push(conf);

					if (printIntermediateConfsToConsole) {
						System.out.println();
						consolePrinter.print(conf, confEcalc.confSpace);
					}
					logPrinter.print(conf, confEcalc.confSpace);

				} else {
					break;
				}
			}

			System.out.println(String.format("Found %d total conformations in energy window", confs.size()));
		}

		return confs;
	}
}
