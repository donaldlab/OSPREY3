package edu.duke.cs.osprey.coffee;

import com.hazelcast.replicatedmap.ReplicatedMap;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.Progress;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;


/**
 * A boltzmann-weighted energy matrix that can be computed in parallel on a cluster.
 */
public class ClusterZMatrix {

	public final ConfEnergyCalculator ecalc;
	public final PosInterGen posInterGen;
	public final BoltzmannCalculator bcalc;


	private final ZMat zmat;

	private class ZMat extends TupleMatrixGeneric<BigExp> {

		BigExp staticStatic;

		private ZMat() {
			super(ecalc.confSpace());
		}

		int numSingles() {
			return getNumOneBody();
		}

		int numPairs() {
			return getNumPairwise();
		}
	}

	public ClusterZMatrix(ConfEnergyCalculator ecalc, PosInterGen posInterGen, BoltzmannCalculator bcalc) {

		this.ecalc = ecalc;
		this.posInterGen = posInterGen;
		this.bcalc = bcalc;

		zmat = new ZMat();
	}

	private BigExp z(double energy) {
		return new BigExp(bcalc.calcPrecise(energy));
	}

	public void compute(ClusterMember member, TaskExecutor tasks) {

		computeStaticStatic(member, tasks);
		computeSingles(member, tasks);
		computePairs(member, tasks);
	}

	private void computeStaticStatic(ClusterMember member, TaskExecutor tasks) {

		// make a replicated map, so hazelcast will sync everything to all members for us
		member.log0("computing static-static ...");
		member.barrier(1, TimeUnit.MINUTES);

		ReplicatedMap<Integer,BigExp> map = member.inst.getReplicatedMap("zmat");

		// only one energy to compute, so do it on member 0
		if (member.id() == 0) {
			Batch batch = new Batch();
			batch.addStaticStatic(0);
			batch.submit(tasks, map, null);
			tasks.waitForFinish();
		}

		// wait for all the members to finish
		member.log0("synchronizing nodes ...");
		member.barrier(5, TimeUnit.MINUTES);

		// copy the replicated map into the local zmat
		zmat.staticStatic = getOrWait(member, map, 0);

		member.barrier(1, TimeUnit.MINUTES);

		// cleanup
		map.destroy();

		member.log0("static-static finished");
	}
	
	private void computeSingles(ClusterMember member, TaskExecutor tasks) {

		// make a replicated map, so hazelcast will sync everything to all members for us
		member.log0("computing %d singles ...", zmat.numSingles());
		member.barrier(1, TimeUnit.MINUTES);

		ReplicatedMap<Integer,BigExp> map = member.inst.getReplicatedMap("zmat");

		var range = member.simplePartition(zmat.numSingles());

		// track progress on the 0 member
		Progress progress = null;
		if (member.id() == 0) {
			progress = new Progress(range.size());
		}

		// statically partition the workload among the members
		Batch batch = new Batch();
		int index = 0;
		for (int posi=0; posi<ecalc.confSpace().numPos(); posi++) {
			for (int confi=0; confi<ecalc.confSpace().numConf(posi); confi++) {
				if (range.contains(index)) {
					batch.addSingle(index, posi, confi);
					if (batch.isFull()) {
						batch.submit(tasks, map, progress);
						batch = new Batch();
					}
				}
				index += 1;
			}
		}
		if (!batch.isEmpty()) {
			batch.submit(tasks, map, progress);
		}
		tasks.waitForFinish();

		// wait for all the members to finish
		member.log0("synchronizing nodes ...");
		member.barrier(5, TimeUnit.MINUTES);

		// copy the replicated map into the local zmat
		index = 0;
		for (int posi=0; posi<ecalc.confSpace().numPos(); posi++) {
			for (int confi=0; confi<ecalc.confSpace().numConf(posi); confi++) {
				zmat.setOneBody(posi, confi, getOrWait(member, map, index));
				index += 1;
			}
		}

		member.barrier(1, TimeUnit.MINUTES);

		// cleanup
		map.destroy();

		member.log0("singles finished");
	}

	private void computePairs(ClusterMember member, TaskExecutor tasks) {

		// make a replicated map, so hazelcast will sync everything to all members for us
		member.log0("computing %d pairs ...", zmat.numPairs());
		member.barrier(1, TimeUnit.MINUTES);
		ReplicatedMap<Integer,BigExp> map = member.inst.getReplicatedMap("zmat");

		var range = member.simplePartition(zmat.numPairs());

		// track progress on the 0 member
		Progress progress = null;
		if (member.id() == 0) {
			progress = new Progress(range.size());
		}

		// statically partition the workload among the members
		Batch batch = new Batch();
		int index = 0;
		for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				for (int confi1=0; confi1<ecalc.confSpace().numConf(posi1); confi1++) {
					for (int confi2=0; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
						if (range.contains(index)) {
							batch.addPair(index, posi1, confi1, posi2, confi2);
							if (batch.isFull()) {
								batch.submit(tasks, map, progress);
								batch = new Batch();
							}
						}
						index += 1;
					}
				}
			}
		}
		if (!batch.isEmpty()) {
			batch.submit(tasks, map, progress);
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
						zmat.setPairwise(posi1, confi1, posi2, confi2, getOrWait(member, map, index));
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
			ClusterMember.sleep(500);
			z = map.get(index);
		}

		return z;
	}

	private class Batch {

		final int capacity;
		final List<Integer> indices;
		final List<ConfEnergyCalculator.MinimizationJob> jobs;

		Batch() {
			capacity = ecalc.maxBatchSize();
			indices = new ArrayList<>(capacity);
			jobs = new ArrayList<>(capacity);
		}

		void addStaticStatic(int index) {
			indices.add(index);
			jobs.add(new ConfEnergyCalculator.MinimizationJob(
				ecalc.confSpace().assign(),
				posInterGen.staticStatic()
			));
		}

		void addSingle(int index, int posi, int confi) {
			indices.add(index);
			jobs.add(new ConfEnergyCalculator.MinimizationJob(
				ecalc.confSpace().assign(posi, confi),
				posInterGen.single(ecalc.confSpace(), posi, confi)
			));
		}

		void addPair(int index, int posi1, int confi1, int posi2, int confi2) {
			indices.add(index);
			jobs.add(new ConfEnergyCalculator.MinimizationJob(
				ecalc.confSpace().assign(posi1, confi1, posi2, confi2),
				posInterGen.pair(ecalc.confSpace(), posi1, confi1, posi2, confi2)
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

		void submit(TaskExecutor tasks, Map<Integer,BigExp> map, Progress progress) {
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
		return zmat.staticStatic;
	}

	public BigExp single(int posi, int confi) {
		return zmat.getOneBody(posi, confi);
	}

	public BigExp pair(int posi1, int confi1, int posi2, int confi2) {
		return zmat.getPairwise(posi1, confi1, posi2, confi2);
	}

	public BigExp pairUpper(int posi1, int confi1, int posi2) {
		// TODO: cache these values
		BigExp max = null;
		for (int confi2=1; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
			BigExp z = zmat.getPairwise(posi1, confi1, posi2, confi2);
			if (max == null || z.greaterThan(max)) {
				max = z;
			}
		}
		return max;
	}
}
