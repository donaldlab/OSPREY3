package edu.duke.cs.osprey.coffee;

import com.hazelcast.replicatedmap.ReplicatedMap;
import edu.duke.cs.osprey.confspace.TripleMatrix;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Progress;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;


/**
 * A boltzmann-weighted energy matrix that can be computed in parallel on a cluster.
 */
public class ClusterZMatrix {

	public final ConfSpace confSpace;
	public final PosInterGen posInterGen;
	public final BoltzmannCalculator bcalc;

	private BigExp staticStatic;
	private final TupleMatrixGeneric<BigExp> singlesPairs;
	private TripleMatrix<BigExp> triples;

	private final BigExp[][] optimizationCache;

	public ClusterZMatrix(ConfSpace confSpace, PosInterGen posInterGen, BoltzmannCalculator bcalc) {

		this.confSpace = confSpace;
		this.posInterGen = posInterGen;
		this.bcalc = bcalc;

		staticStatic = new BigExp(1.0, 0);
		singlesPairs = new TupleMatrixGeneric<>(confSpace);
		triples = null;

		// allocate memory for the optimization cache
		optimizationCache = new BigExp[singlesPairs.getNumPairwise()][];
		for (int posi1=0; posi1<confSpace.numPos(); posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				optimizationCache[singlesPairs.getPairwiseIndex(posi1, posi2)] = new BigExp[confSpace.numConf(posi1)];
			}
		}
	}

	private BigExp z(double energy) {
		return new BigExp(bcalc.calcPrecise(energy));
	}

	public void compute(ClusterMember member, TaskExecutor tasks, boolean includeStaticStatic, Double tripleCorrectionThreshold, ConfEnergyCalculator ecalc) {

		if (includeStaticStatic) {
			computeStaticStatic(member, tasks, ecalc);
		}
		computeSingles(member, tasks, ecalc);
		computePairs(member, tasks, ecalc);
		if (tripleCorrectionThreshold != null) {
			computeTripleCorrections(member, tasks, ecalc, tripleCorrectionThreshold);
		}
	}

	private void computeStaticStatic(ClusterMember member, TaskExecutor tasks, ConfEnergyCalculator ecalc) {

		// make a replicated map, so hazelcast will sync everything to all members for us
		member.log0("computing static-static ...");
		member.barrier(1, TimeUnit.MINUTES);

		ReplicatedMap<Integer,BigExp> map = member.inst.getReplicatedMap("zmat");

		// only one energy to compute, so do it on member 0
		if (member.id() == 0) {
			Batch batch = new Batch(ecalc.maxBatchSize());
			batch.addStaticStatic(0);
			batch.submit(tasks, map, null, ecalc);
			tasks.waitForFinish();
		}

		// wait for all the members to finish
		member.log0("synchronizing nodes ...");
		member.barrier(5, TimeUnit.MINUTES);

		// copy the replicated map into the local zmat
		staticStatic = getOrWait(member, map, 0);

		member.barrier(1, TimeUnit.MINUTES);

		// cleanup
		map.destroy();

		member.log0("static-static finished");
	}

	private void computeSingles(ClusterMember member, TaskExecutor tasks, ConfEnergyCalculator ecalc) {

		// make a replicated map, so hazelcast will sync everything to all members for us
		member.log0("computing %d singles ...", singlesPairs.getNumOneBody());
		member.barrier(1, TimeUnit.MINUTES);

		ReplicatedMap<Integer,BigExp> map = member.inst.getReplicatedMap("zmat");

		var range = member.simplePartition(singlesPairs.getNumOneBody());

		// track progress on the 0 member
		Progress progress = null;
		if (member.id() == 0) {
			progress = new Progress(range.size());
		}

		// statically partition the workload among the members
		Batch batch = new Batch(ecalc.maxBatchSize());
		int index = 0;
		for (int posi=0; posi<ecalc.confSpace().numPos(); posi++) {
			for (int confi=0; confi<ecalc.confSpace().numConf(posi); confi++) {
				if (range.contains(index)) {
					batch.addSingle(index, posi, confi);
					if (batch.isFull()) {
						batch.submit(tasks, map, progress, ecalc);
						batch = new Batch(ecalc.maxBatchSize());
					}
				}
				index += 1;
			}
		}
		if (!batch.isEmpty()) {
			batch.submit(tasks, map, progress, ecalc);
		}
		tasks.waitForFinish();

		// wait for all the members to finish
		member.log0("synchronizing nodes ...");
		member.barrier(5, TimeUnit.MINUTES);

		// copy the replicated map into the local zmat
		index = 0;
		for (int posi=0; posi<ecalc.confSpace().numPos(); posi++) {
			for (int confi=0; confi<ecalc.confSpace().numConf(posi); confi++) {
				singlesPairs.setOneBody(posi, confi, getOrWait(member, map, index));
				index += 1;
			}
		}

		member.barrier(1, TimeUnit.MINUTES);

		// cleanup
		map.destroy();

		member.log0("singles finished");
	}

	private void computePairs(ClusterMember member, TaskExecutor tasks, ConfEnergyCalculator ecalc) {

		// make a replicated map, so hazelcast will sync everything to all members for us
		member.log0("computing %d pairs ...", singlesPairs.getNumPairwise());
		member.barrier(1, TimeUnit.MINUTES);
		ReplicatedMap<Integer,BigExp> map = member.inst.getReplicatedMap("zmat");

		var range = member.simplePartition(singlesPairs.getNumPairwise());

		// track progress on the 0 member
		Progress progress = null;
		if (member.id() == 0) {
			progress = new Progress(range.size());
		}

		// statically partition the workload among the members
		Batch batch = new Batch(ecalc.maxBatchSize());
		int index = 0;
		for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				for (int confi1=0; confi1<ecalc.confSpace().numConf(posi1); confi1++) {
					for (int confi2=0; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
						if (range.contains(index)) {
							batch.addPair(index, posi1, confi1, posi2, confi2);
							if (batch.isFull()) {
								batch.submit(tasks, map, progress, ecalc);
								batch = new Batch(ecalc.maxBatchSize());
							}
						}
						index += 1;
					}
				}
			}
		}
		if (!batch.isEmpty()) {
			batch.submit(tasks, map, progress, ecalc);
		}
		tasks.waitForFinish();

		// wait for all the members to finish
		member.log0("synchronizing nodes ...");
		member.barrier(1, TimeUnit.MINUTES);

		// copy the replicated map into the local zmat
		index = 0;
		for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				for (int confi1=0; confi1<ecalc.confSpace().numConf(posi1); confi1++) {
					for (int confi2=0; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
						singlesPairs.setPairwise(posi1, confi1, posi2, confi2, getOrWait(member, map, index));
						index += 1;
					}
				}
			}
		}

		member.barrier(1, TimeUnit.MINUTES);

		// cleanup
		map.destroy();

		member.log0("pairs finished");
	}

	private void computeTripleCorrections(ClusterMember member, TaskExecutor tasks, ConfEnergyCalculator ecalc, double energyThreshold) {

		// allocate space
		triples = new TripleMatrix<>(confSpace);

		// filter the triples by the threshold
		int numTriples = 0;
		BigExp zThreshold = new BigExp(bcalc.calcPrecise(energyThreshold));
		boolean[] passed = new boolean[triples.size()];
		Arrays.fill(passed, false);
		int index = 0;
		for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
			int n1 = ecalc.confSpace().numConf(posi1);
			for (int posi2=0; posi2<posi1; posi2++) {
				int n2 = ecalc.confSpace().numConf(posi2);
				for (int posi3=0; posi3<posi2; posi3++) {
					int n3 = ecalc.confSpace().numConf(posi3);

					for (int confi1=0; confi1<n1; confi1++) {

						// skip triples whose constituents are below the z threshold
						if (singlesPairs.getOneBody(posi1, confi1).lessThan(zThreshold)) {
							index += n2*n3;
							continue;
						}

						for (int confi2=0; confi2<n2; confi2++) {

							// skip triples whose constituents are below the z threshold
							if (singlesPairs.getOneBody(posi2, confi2).lessThan(zThreshold)
							|| singlesPairs.getPairwise(posi1, confi1, posi2, confi2).lessThan(zThreshold)) {
								index += n3;
								continue;
							}

							for (int confi3=0; confi3<n3; confi3++) {

								// skip triples whose constituents are below the z threshold
								if (singlesPairs.getOneBody(posi3, confi3).lessThan(zThreshold)
								|| singlesPairs.getPairwise(posi1, confi1, posi3, confi3).lessThan(zThreshold)
								|| singlesPairs.getPairwise(posi2, confi2, posi3, confi3).lessThan(zThreshold)) {
									index += 1;
									continue;
								}

								passed[index++] = true;
								numTriples += 1;
							}
						}
					}
				}
			}
		}

		// make a replicated map, so hazelcast will sync everything to all members for us
		member.log0("computing %d triples ...", numTriples);
		member.barrier(1, TimeUnit.MINUTES);
		ReplicatedMap<Integer,BigExp> map = member.inst.getReplicatedMap("zmat");

		var range = member.simplePartition(numTriples);

		// track progress on the 0 member
		Progress progress = null;
		if (member.id() == 0) {
			progress = new Progress(range.size());
		}

		// statically partition the workload among the members
		Batch batch = new Batch(ecalc.maxBatchSize());
		index = 0;
		int triplei = 0;
		for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				for (int posi3=0; posi3<posi2; posi3++) {
					for (int confi1=0; confi1<ecalc.confSpace().numConf(posi1); confi1++) {
						for (int confi2=0; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
							for (int confi3=0; confi3<ecalc.confSpace().numConf(posi3); confi3++) {
								if (passed[index]) {
									if (range.contains(triplei)) {
										batch.addTripleCorrection(index, posi1, confi1, posi2, confi2, posi3, confi3);
										if (batch.isFull()) {
											batch.submit(tasks, map, progress, ecalc);
											batch = new Batch(ecalc.maxBatchSize());
										}
										triplei += 1;
									}
								}
								index += 1;
							}
						}
					}
				}
			}
		}
		if (!batch.isEmpty()) {
			batch.submit(tasks, map, progress, ecalc);
		}
		tasks.waitForFinish();

		// wait for all the members to finish
		member.log0("synchronizing nodes ...");
		member.barrier(1, TimeUnit.MINUTES);

		// copy the replicated map into the local zmat
		index = 0;
		for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				for (int posi3=0; posi3<posi2; posi3++) {
					for (int confi1=0; confi1<ecalc.confSpace().numConf(posi1); confi1++) {
						for (int confi2=0; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
							for (int confi3=0; confi3<ecalc.confSpace().numConf(posi3); confi3++) {

								if (passed[index]) {

									BigExp triple = getOrWait(member, map, index);

									// convert the triple energy into a correction
									BigExp divisor = new BigExp(pair(posi1, confi1, posi2, confi2));
									divisor.mult(pair(posi1, confi1, posi3, confi3));
									divisor.mult(pair(posi2, confi2, posi3, confi3));
									divisor.pow(1.0/MathTools.numTriplesPerPair(ecalc.confSpace().numPos()));

									// only use unhelpful corrections
									if (triple.lessThan(divisor)) {
										triple.div(divisor);
										triples.set(posi1, confi1, posi2, confi2, posi3, confi3, triple);
									}
								}
								index += 1;
							}
						}
					}
				}
			}
		}

		member.barrier(1, TimeUnit.MINUTES);

		// cleanup
		map.destroy();

		member.log0("triples finished");
	}

	private BigExp getOrWait(ClusterMember member, Map<Integer,BigExp> map, Integer index) {

		BigExp z = map.get(index);

		long maxWaitNs = 60*1_000_000_000L; // 1 min
		long startNs = System.nanoTime();

		while (z == null) {

			// value not ready yet, wait a bit and try again
			long elapsedNs = System.nanoTime() - startNs;
			if (elapsedNs > maxWaitNs) {
				member.throwTimeout("timed out waiting for matrix element from cluster");
			}
			ThreadTools.sleep(500);
			z = map.get(index);
		}

		return z;
	}

	private class Batch {

		final int capacity;
		final List<Integer> indices;
		final List<ConfEnergyCalculator.MinimizationJob> jobs;

		Batch(int capacity) {
			this.capacity = capacity;
			indices = new ArrayList<>(capacity);
			jobs = new ArrayList<>(capacity);
		}

		void addStaticStatic(int index) {
			indices.add(index);
			jobs.add(new ConfEnergyCalculator.MinimizationJob(
				confSpace.assign(),
				posInterGen.staticStatic()
			));
		}

		void addSingle(int index, int posi, int confi) {
			indices.add(index);
			jobs.add(new ConfEnergyCalculator.MinimizationJob(
				confSpace.assign(posi, confi),
				posInterGen.single(confSpace, posi, confi)
			));
		}

		void addPair(int index, int posi1, int confi1, int posi2, int confi2) {
			indices.add(index);
			jobs.add(new ConfEnergyCalculator.MinimizationJob(
				confSpace.assign(posi1, confi1, posi2, confi2),
				posInterGen.pair(confSpace, posi1, confi1, posi2, confi2)
			));
		}

		void addTripleCorrection(int index, int posi1, int confi1, int posi2, int confi2, int posi3, int confi3) {
			indices.add(index);
			jobs.add(new ConfEnergyCalculator.MinimizationJob(
				confSpace.assign(posi1, confi1, posi2, confi2, posi3, confi3),
				posInterGen.tripleCorrection(confSpace, posi1, confi1, posi2, confi2, posi3, confi3)
			));
		}

		int size() {
			return jobs.size();
		}

		public boolean isEmpty() {
			return size() == 0;
		}

		boolean isFull() {
			return size() == capacity;
		}

		void submit(TaskExecutor tasks, Map<Integer,BigExp> map, Progress progress, ConfEnergyCalculator ecalc) {
			assert (!isEmpty());
			tasks.submit(
				() -> {
					ecalc.minimizeEnergies(jobs);
					return jobs.stream()
						.map(job -> z(job.energy))
						.collect(Collectors.toList());
				},
				zs -> {
					for (int j=0; j<size(); j++) {
						map.put(indices.get(j), zs.get(j));
					}
					if (progress != null) {
						assert (size() != 0);
						progress.incrementProgress(size());
					}
				}
			);
		}
	}

	public BigExp staticStatic() {
		return staticStatic;
	}

	public BigExp single(int posi, int confi) {
		return singlesPairs.getOneBody(posi, confi);
	}

	public BigExp pair(int posi1, int confi1, int posi2, int confi2) {
		return singlesPairs.getPairwise(posi1, confi1, posi2, confi2);
	}

	public BigExp pairUpper(int posi1, int confi1, int posi2) {

		// check the cache first
		int i = singlesPairs.getPairwiseIndex(posi1, posi2);
		BigExp max = optimizationCache[i][confi1];
		if (max != null) {
			return max;
		}

		// cache miss, do the optimization
		for (int confi2=1; confi2<confSpace.numConf(posi2); confi2++) {
			BigExp z = singlesPairs.getPairwise(posi1, confi1, posi2, confi2);
			if (max == null || z.greaterThan(max)) {
				max = z;
			}
		}

		optimizationCache[i][confi1] = max;

		return max;
	}

	public boolean hasTriples() {
		return triples != null;
	}

	public BigExp triple(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3) {
		return triples.get(posi1, confi1, posi2, confi2, posi3, confi3);
	}
}
