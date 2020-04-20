package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.HigherOrderGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.JvmMem;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.Streams;

import java.io.File;
import java.io.Serializable;
import java.time.LocalDateTime;
import java.util.*;
import java.util.concurrent.CountDownLatch;
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
	public void doEeet(Sequence seq, int numConfs, int outputIntervalSeconds, String pdbPath) {

		long outputIntervalNs = outputIntervalSeconds * 1_000_000_000L;

		try (TaskExecutor.ContextGroup contexts = tasks.contextGroup()) {

			// set up task contexts
			contexts.putContext(0, EnergyTask.class, new TaskContext(confEcalc));
			contexts.putContext(0, ConfTask.class, new TaskContext(confEcalc));

			// skip the calculation on member nodes
			if (tasks instanceof Cluster.Member) {
				return;
			}

			// make A* go BRRRRRRR
			RCs rcs = seq.makeRCs(confEcalc.confSpace);
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
				.setCustom(
					new DynamicHMeanAStarOrder(MathTools.Optimizer.Minimize),
					new HigherOrderGScorer(emat),
					// TODO: the H-scorer is a bottleneck, needs to be optimized
					//new TriplewiseHScorer(emat)
					// luckily, the usual pairwise H-scorer is still a decently fast heuristic,
					// even with the higher-order G-scorer!
					// but it will cause A* to eat more memory than usual =(
					new TraditionalPairwiseHScorer(emat, rcs)
				)
				.build();

			Stats stats = new Stats();

			log("Searching for %d lowest-energy conformations", numConfs);

			Stopwatch stopwatch = new Stopwatch().start();

			try (ThreadPoolTaskExecutor ensembleTasks = new ThreadPoolTaskExecutor()) {
				ensembleTasks.start(1);

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
							if (stats.needsOutput && nowNs > stats.lastOutputNs + outputIntervalNs) {

								// copy the confs before giving them to another thread
								var bestConfs = new ArrayList<>(stats.bestConfs.values());

								ensembleTasks.submit(
									() -> {
										writeEnsemble(seq, bestConfs, pdbPath);
										return 42; // it's the answer
									},
									answer -> {
										synchronized (stats) { // don't race the main thread
											stats.lastOutputNs = System.nanoTime();
											stats.needsOutput = true;
										}
									}
								);
								stats.needsOutput = false;
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
								while (stats.bestConfs.size() > numConfs) {
									stats.bestConfs.pollLastEntry();
								}

								// show every minimization?
								if (false) {
									log("conf %10d   score %12.6f   energy %12.6f   gap %12.6f",
										stats.numConfsEnergied,
										econf.getScore(),
										econf.getEnergy(),
										econf.getEnergy() - econf.getScore()
									);
								}

								// should we show progress?
								long nowNs = System.nanoTime();
								final double reportIntervalS = 5.0;
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
								boolean isDone = stats.bestConfs.size() >= numConfs
									&& stats.bestConfs.keySet().stream().allMatch(energy -> energy < stats.maxScore);
								if (isDone) {
									keepSearching.set(false);
								}
							}
						}
					);
				}

				tasks.waitForFinish();
				ensembleTasks.waitForFinish();

				log("Found %d lowest-energy confs in %s! We're done here!",
					numConfs, stopwatch.stop().getTime(2)
				);

				// write the ensemble one final time
				writeEnsemble(seq, stats.bestConfs.values(), pdbPath);
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
		boolean needsOutput = true;
	}

	private static class TaskContext {

		final ConfEnergyCalculator confEcalc;

		public TaskContext(ConfEnergyCalculator confEcalc) {
			this.confEcalc = confEcalc;
		}
	}

	private static class EnergyTask extends Cluster.Task<ConfSearch.EnergiedConf,TaskContext> {

		final ConfSearch.ScoredConf conf;

		EnergyTask(ConfSearch.ScoredConf conf) {
			super(0);
			this.conf = conf;
		}

		@Override
		public ConfSearch.EnergiedConf run(TaskContext ctx) {
			return ctx.confEcalc.calcEnergy(conf);
		}
	}

	private static class ConfResult implements Serializable {

		final double[][] coords;
		final double energy;

		public ConfResult(EnergyCalculator.EnergiedParametricMolecule epmol) {

			// gather all the residue coords, for simple serialization
			coords = new double[epmol.pmol.mol.residues.size()][];
			for (Residue res : epmol.pmol.mol.residues) {
				coords[res.indexInMolecule] = res.coords;
			}

			this.energy = epmol.energy;
		}

		public EnergyCalculator.EnergiedParametricMolecule makeEpmol(ParametricMolecule pmol) {
			for (Residue res : pmol.mol.residues) {
				res.coords = coords[res.indexInMolecule];
			}
			return new EnergyCalculator.EnergiedParametricMolecule(pmol, null, energy);
		}
	}

	private static class ConfTask extends Cluster.Task<ConfResult,TaskContext> {

		final ConfSearch.ScoredConf conf;

		ConfTask(ConfSearch.ScoredConf conf) {
			super(0);
			this.conf = conf;
		}

		@Override
		public ConfResult run(TaskContext ctx) {
			return new ConfResult(ctx.confEcalc.calcEnergy(new RCTuple(conf.getAssignments())));
		}
	}

	private void writeEnsemble(Sequence seq, Collection<ConfSearch.EnergiedConf> econfs, String pdbPath) {

		log("Collecting latest ensemble ...");

		CountDownLatch latch = new CountDownLatch(econfs.size());

		// minimize all the conformations (again) to get the coords
		List<EnergyCalculator.EnergiedParametricMolecule> epmols = new ArrayList<>(econfs.size());
		for (ConfSearch.EnergiedConf econf : econfs) {
			tasks.submit(
				new ConfTask(econf),
				confResult -> {
					epmols.add(confResult.makeEpmol(confEcalc.confSpace.makeMolecule(econf)));
					latch.countDown();
				}
			);
		}

		// wait for the minimizations to finish
		// NOTE: can't use tasks.waitForFinish() here, since the main thread is still using tasks
		try {
			latch.await();
		} catch (InterruptedException ex) {
			throw new RuntimeException("interrupted waiting for ensemble minimizations", ex);
		}

		// sort the confs by energy
		epmols.sort(Comparator.comparing(epmol -> epmol.energy));

		// get a useful id for the sequence
		String seqId = Streams.of(seq.assignments())
			.filter(assignment -> assignment.isMutated())
			.map(assignment -> ResNum.render(assignment) + "=" + ResTypeMutations.render(assignment))
			.collect(Collectors.joining(" "));
		if (seqId.isEmpty()) {
			seqId = "wild-type";
		}

		// write the PDB
        LocalDateTime currTime = LocalDateTime.now();
		String nowDirName = String.format("%02d-%02d-%02d", currTime.getHour(), currTime.getMinute(), currTime.getSecond());
		File nowDir = new File(String.format("%s/%s", pdbPath, nowDirName));
		if(!nowDir.exists())
			nowDir.mkdir();
		File file = new File(String.format("%s/%s/seq.%s.pdb", pdbPath, nowDirName, seqId));
		PDBIO.writeFile(
			epmols, file,
			String.format("Ensemble of %d lowest-energy conformations for sequence: %s",
				epmols.size(), seq.toString(Sequence.Renderer.AssignmentMutations)
			)
		);

		log("Wrote sequences ensemble to %s", file.getAbsolutePath());
	}
}
