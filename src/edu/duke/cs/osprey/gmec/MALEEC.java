package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.JvmMem;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.Streams;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.confspace.Sequence.Renderer.ResNum;
import static edu.duke.cs.osprey.confspace.Sequence.Renderer.ResTypeMutations;
import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Massively Awesome Low-Energy Ensemble Calculator
 */
public class MALEEC {

	public final TaskExecutor tasks;
	public final ConfEnergyCalculator confEcalc;
	public final EnergyMatrix emat;

	public MALEEC(TaskExecutor tasks, ConfEnergyCalculator confEcalc, EnergyMatrix emat) {
		this.tasks = tasks;
		this.confEcalc = confEcalc;
		this.emat = emat;
	}

	/**
	 * Start a GMEC-like search to find the lowest N conformations by minimized energy.
	 *
	 * Every S seconds, write out an ensemble of the best conformations known so far.
	 *
	 * LEEEEEEROOYYYY!!!
	 */
	// at least I've got chicken
	public void doEeet(Sequence seq, int maxNumConfs, int outputIntervalSeconds) {

		long outputIntervalNs = outputIntervalSeconds * 1_000_000_000L;

		try (TaskExecutor.ContextGroup contexts = tasks.contextGroup()) {

			// set up task contexts
			contexts.putContext(0, EnergyTask.class, new EnergyTask.Context(confEcalc));

			// skip the calculation on member nodes
			if (tasks instanceof Cluster.Member) {
				return;
			}

			// configure A* to go BRRRRRRR
			RCs rcs = seq.makeRCs(confEcalc.confSpace);
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
				.setTraditional() // TODO: update A* to use triples and quads
				.build();

			Stats stats = new Stats();

			// make an extra thread to collect the ensembles
			try (ThreadPoolTaskExecutor ensembleTasks = new ThreadPoolTaskExecutor()) {
				ensembleTasks.start(1);

				Stopwatch stopwatch = new Stopwatch().start();

				AtomicBoolean keepSearching = new AtomicBoolean(true);
				while (keepSearching.get()) {

					// get the next conf from A*
					ConfSearch.ScoredConf conf = astar.nextConf();
					if (conf == null) {
						log("A* ran out of confs. Hooray!");
						break;
					}

					{
						long nowNs = System.nanoTime();

						// update the stats
						synchronized (stats) { // don't race the listener thread

							stats.numConfsScored += 1;
							stats.maxScore = conf.getScore();

							// write the output ensemble if needed
							if (nowNs > stats.lastOutputNs + outputIntervalNs) {
								writeEnsemble(ensembleTasks, seq, stats.bestConfs.values());
								stats.lastOutputNs = nowNs;
							}
						}
					}

					// calculate the energy
					tasks.submit(
						new EnergyTask(conf),
						econf -> {

							// are we done yet?
							synchronized (stats) { // don't race the main thread

								stats.numConfsEnergied += 1;
								stats.minEnergy = Math.min(stats.minEnergy, econf.getEnergy());

								// update the best confs
								stats.bestConfs.put(econf.getEnergy(), econf);
								while (stats.bestConfs.size() > maxNumConfs) {
									stats.bestConfs.pollLastEntry();
								}

								// should we show progress?
								long nowNs = System.nanoTime();
								final double reportIntervalS = 1.0;
								final long reportIntervalNs = (long)(reportIntervalS * 1e9);
								if (nowNs > stats.lastReportNs + reportIntervalNs) {

									// yup, update the performance logs
									stats.confsEnergiedHistory.add(stats.numConfsEnergied - stats.lastReportConfsEnergied);
									while (stats.confsEnergiedHistory.size() > 30) {
										stats.confsEnergiedHistory.removeFirst();
									}

									double bestEnergy = stats.bestConfs.firstKey();
									double worstEnergy = stats.bestConfs.lastKey();

									// show the progress
									log("best confs: %3d/%10d   %6.1f confs/sec   energy range: [%12.6f,%12.6f]   max score: %12.6f   gap: %12.6f   time:%10s   heapMem:%s",
										stats.bestConfs.size(),
										stats.numConfsEnergied,
										stats.confsEnergiedHistory.stream().mapToLong(i -> i).sum()/reportIntervalS/stats.confsEnergiedHistory.size(),
										bestEnergy,
										worstEnergy,
										stats.maxScore,
										worstEnergy - stats.maxScore,
										stopwatch.getTime(2),
										JvmMem.getOldPool()
									);

									stats.lastReportNs = nowNs;
									stats.lastReportConfsEnergied = stats.numConfsEnergied;
								}

								// if all the best confs have energies below the max score, we're done
								boolean isDone = stats.bestConfs.size() >= maxNumConfs
									&& stats.bestConfs.keySet().stream().allMatch(energy -> energy < stats.maxScore);
								if (isDone) {
									keepSearching.set(false);
								}
							}
						}
					);
				}

				tasks.waitForFinish();

				log("Found %d lowest-energy confs in %s! We're done here!",
					maxNumConfs, stopwatch.stop().getTime(2)
				);

				// write the ensemble one final time
				writeEnsemble(ensembleTasks, seq, stats.bestConfs.values());
				ensembleTasks.waitForFinish();
			}

		} catch (Throwable t) {
			// explicitly show exceptions now, before the ContextGroup close()es and hides them
			t.printStackTrace();
			throw new Error();
		}
	}

	private static class Stats {
		long numConfsScored = 0;
		long numConfsEnergied = 0;
		double maxScore = Double.NEGATIVE_INFINITY;
		double minEnergy = Double.POSITIVE_INFINITY;
		TreeMap<Double,ConfSearch.EnergiedConf> bestConfs = new TreeMap<>();
		long lastReportNs = 0;
		long lastReportConfsEnergied = 0;
		Deque<Long> confsEnergiedHistory = new ArrayDeque<>();
		long lastOutputNs = System.nanoTime();
	}

	private static class EnergyTask extends Cluster.Task<ConfSearch.EnergiedConf,EnergyTask.Context> {

		static class Context {

			final ConfEnergyCalculator confEcalc;

			public Context(ConfEnergyCalculator confEcalc) {
				this.confEcalc = confEcalc;
			}
		}

		final ConfSearch.ScoredConf conf;

		EnergyTask(ConfSearch.ScoredConf conf) {
			super(0);
			this.conf = conf;
		}

		@Override
		public ConfSearch.EnergiedConf run(Context ctx) {
			return ctx.confEcalc.calcEnergy(conf);
		}
	}

	private void writeEnsemble(TaskExecutor tasks, Sequence seq, Collection<ConfSearch.EnergiedConf> econfs) {

		// copy the confs for the task thread
		var econfsCopy = new ArrayList<>(econfs);

		tasks.submit(
			() -> {
				log("Collecting latest ensemble ...");

				// sort the confs by energy
				List<EnergyCalculator.EnergiedParametricMolecule> epmols = econfsCopy.stream()
					.sorted(Comparator.comparing(econf -> econf.getEnergy()))
					.map(econf -> {
						RCTuple tuple = new RCTuple(econf.getAssignments());
						return confEcalc.calcEnergy(tuple);
					})
					.collect(Collectors.toList());

				// get a useful id for the sequence
				String seqId = Streams.of(seq.assignments())
					.filter(assignment -> assignment.isMutated())
					.map(assignment -> ResNum.render(assignment) + "=" + ResTypeMutations.render(assignment))
					.collect(Collectors.joining(" "));
				if (seqId.isEmpty()) {
					seqId = "wild-type";
				}

				// convert the sequence into a filename
				File file = new File(String.format("seq.%s.pdb", seqId));

				// write the PDB
				PDBIO.writeFile(
					epmols, file,
					String.format("Ensemble of %d lowest-energy conformations for sequence: %s",
						econfsCopy.size(), seq.toString(Sequence.Renderer.AssignmentMutations)
					)
				);

				log("Wrote sequences ensemble to %s", file.getAbsolutePath());

				return 42; // it's the answer
			},
			answer -> {} // don't need the listener thread here
		);
	}
}
