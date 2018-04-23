package edu.duke.cs.osprey.lute;


import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
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
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new LUTEGScorer(confEcalc),
				new LUTEHScorer(confEcalc)
			)
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
