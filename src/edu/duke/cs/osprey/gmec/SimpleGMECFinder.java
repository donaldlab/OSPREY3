package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.externalMemory.EnergiedConfFIFOSerializer;
import edu.duke.cs.osprey.externalMemory.EnergiedConfPrioritySerializer;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.externalMemory.ScoredConfFIFOSerializer;
import edu.duke.cs.osprey.gmec.GMECFinder.ConfPruner;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.Progress;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.*;
import java.util.*;
import java.util.function.Consumer;

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
		 * If a design experiences an unexpected abort, the resume log can allow you to restore the
		 * design state and resume the calculation close to where it was aborted. Set a file to turn on the log.
		 */
		protected File resumeLog = null;
		
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

		public Builder setResumeLog(File val) {
			this.resumeLog = val;
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
				resumeLog
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
	private final ResumeLog resumeLog;
	
	protected SimpleGMECFinder(ConfSearch search, ConfEnergyCalculator confEcalc, ConfPruner pruner, ConfPrinter logPrinter, ConfPrinter consolePrinter, boolean printIntermediateConfsToConsole, boolean useExternalMemory, File resumeLogFile) {
		this.search = search;
		this.confEcalc = confEcalc;
		this.pruner = pruner;
		this.logPrinter = logPrinter;
		this.consolePrinter = consolePrinter;
		this.printIntermediateConfsToConsole = printIntermediateConfsToConsole;
		
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

		this.resumeLog = new ResumeLog(resumeLogFile);
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
		
		// evaluate the min score conf
		System.out.println("Computing energy of min score conf...");
		EnergiedConf eMinScoreConf = confEcalc.calcEnergy(minScoreConf);
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
		
			checkMoreConfs(unpeekedConfs, erange, econfs);
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
	}
	
	private void checkMoreConfs(ConfSearch search, EnergyRange erange, Queue<EnergiedConf> econfs) {
		
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
				EnergiedConf econf = confEcalc.calcEnergy(otherLowEnergyConfs.poll());
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
			
			minimizeLowEnergyConfs(otherLowEnergyConfs, erange, econfs);
		}
	}

	private void minimizeLowEnergyConfs(Queue.FIFO<ScoredConf> lowEnergyConfs, EnergyRange erange, Queue<EnergiedConf> econfs) {

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

		// do we have a resume log?
		if (resumeLog.hasLog()) {

			System.out.println("Restoring state from resume log: " + resumeLog.logFile);

			// yup, read in the already-computed energies
			resumeLog.readAll(lowEnergyConfs, (econf) -> {
				ecalcListener.onFinished(econf);
			});

			System.out.println(String.format("Restored %d conformations. Resuming design...", econfs.size()));

			// TODO: reset progress history, so we get an accurate ETA on the remaining work

		} else if (resumeLog.isLogging()) {

			System.out.println("Starting resume log: " + resumeLog.logFile);
		}

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

			// update the resume log if needed
			if (resumeLog.isLogging()) {
				synchronized (resumeLog) {
					resumeLog.expect(conf);
				}
			}

			// send the conf to the energy calculator
			confEcalc.calcEnergyAsync(conf, (econf) -> {

				// NOTE: this is called on a listener thread, which is separate from the main thread

				// update the resume log if needed
				if (resumeLog.isLogging()) {
					synchronized (resumeLog) {
						resumeLog.save(econf);
					}
				}

				ecalcListener.onFinished(econf);
			});
		}

		confEcalc.tasks.waitForFinish();

		// finished! cleanup the log
		if (resumeLog.isLogging()) {
			System.out.println("GMEC search complete, cleaning up resume log: " + resumeLog.logFile);
			resumeLog.delete();
		}
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

	private class ResumeLog {

		private static final int FlushEveryNConfs = 8;

		public final File logFile;

		private Deque<ScoredConf> expectedConfs = new ArrayDeque<>();
		private List<Double> energies = new ArrayList<>();
		private IntEncoding encoding = null;

		public ResumeLog(File logFile) {

			this.logFile = logFile;

			if (isLogging()) {

				// determine the best encoding scheme for RCs
				encoding = IntEncoding.get(confEcalc.confSpace.positions.stream()
					.map((pos) -> pos.resConfs.size())
					.max(Integer::compare)
					.orElseThrow(() -> new IllegalStateException("conf space has no positions"))
				);
			}
		}

		public boolean isLogging() {
			return logFile != null;
		}

		public boolean hasLog() {
			return isLogging() && logFile.exists();
		}

		public void expect(ScoredConf conf) {
			expectedConfs.add(conf);
		}

		private int findExpectedIndex(ScoredConf conf) {

			// linear search should be plenty fast enough here
			// the expected list should always be pretty small
			int i = 0;
			for (ScoredConf expectedConf : expectedConfs) {
				if (Arrays.equals(conf.getAssignments(), expectedConf.getAssignments())) {
					return i;
				}
				i++;
			}

			throw new NoSuchElementException("expected conf " + Arrays.toString(conf.getAssignments()) + " not found");
		}

		public void save(EnergiedConf econf) {

			int i = findExpectedIndex(econf);

			// expand the energy buffer to make space
			while (energies.size() <= i) {
				energies.add(null);
			}

			// save the energy in the buffer
			energies.set(i, econf.getEnergy());

			// save energies to file as needed
			if (isEnergyBufferFull()) {
				flushEnergies();
			}
		}

		private boolean isEnergyBufferFull() {

			// do we have enough expected confs?
			if (expectedConfs.size() < FlushEveryNConfs) {
				return false;
			}

			// do we have enough consecutive energies?
			if (energies.size() < FlushEveryNConfs) {
				return false;
			}
			for (int i=0; i<FlushEveryNConfs; i++) {
				if (energies.get(i) == null) {
					return false;
				}
			}

			return true;
		}

		private void flushEnergies() {

			try (DataOutputStream out = new DataOutputStream(new FileOutputStream(logFile, true))) {

				for (int i=0; i<energies.size(); i++) {

					// get the next energy to save
					Double energy = energies.get(i);
					if (energy == null) {
						break;
					}

					// get the corresponding conf
					ScoredConf conf = expectedConfs.pop();

					// save it to the file
					for (int rc : conf.getAssignments()) {
						encoding.write(out, rc);
					}
					out.writeDouble(conf.getScore());
					out.writeDouble(energy);
				}

			} catch (IOException ex) {
				throw new RuntimeException("Can't write to resume log: " + logFile, ex);
			}
		}

		public void readAll(Queue<ScoredConf> confs, Consumer<EnergiedConf> onEconf) {

			List<SimpleConfSpace.Position> positions = confEcalc.confSpace.positions;

			try (FileInputStream fin = new FileInputStream(logFile)) {
				DataInputStream in = new DataInputStream(fin);

				while (true) {

					// are we at the end of the file?
					if (fin.getChannel().position()  == fin.getChannel().size()) {

						// end of the log, we're done reading
						break;
					}

					// read the next log conf
					int[] assignments = new int[positions.size()];
					for (SimpleConfSpace.Position pos : positions) {
						assignments[pos.index] = encoding.read(in);
					}
					double score = in.readDouble();
					double energy = in.readDouble();
					EnergiedConf econf = new EnergiedConf(assignments, score, energy);

					// compare to the next expected conf
					if (confs.isEmpty()) {
						throw new IllegalStateException(String.format("resume log expected conformation %f %s, but didn't find any.",
							econf.getScore(), Arrays.toString(econf.getAssignments())
						));
					}
					ScoredConf conf = confs.poll();

					if (econf.getScore() != conf.getScore() || !Arrays.equals(econf.getAssignments(), conf.getAssignments())) {
						throw new IllegalStateException(String.format("resume log expected conformation %f %s, but but found %f %s instead.",
							econf.getScore(), Arrays.toString(econf.getAssignments()),
							conf.getScore(), Arrays.toString(conf.getAssignments())
						));
					}

					// all is well, pass up the logged conf
					onEconf.accept(econf);
				}
			} catch (IOException ex) {
				throw new RuntimeException("Can't read from resume log: " + logFile, ex);
			}
		}

		public void delete() {
			logFile.delete();
		}
	}
}
