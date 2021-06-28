package edu.duke.cs.osprey.coffee.zmat;

import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.confspace.TripleMatrix;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Generator;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.WorkLatch;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.Progress;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;


/**
 * A boltzmann-weighted energy matrix that can be computed in parallel on a cluster.
 */
public class ClusterZMatrix {

	public static final String ServiceName = "ClusterZMatrix";

	public final ConfSpace confSpace;
	public final PosInterGen posInterGen;
	public final BoltzmannCalculator bcalc;

	private BigExp staticStatic;
	private final TupleMatrixGeneric<BigExp> singlesPairs;
	private TripleMatrix<BigExp> triples;

	private WorkLatch latch = null;
	private Progress progress = null;
	private AtomicLong droppedTuples = null;

	// Allow this class to compute a "rigid zmatrix" used for computing pfunc lower bounds
	boolean lowerBounder;

	public ClusterZMatrix(ConfSpace confSpace, PosInterGen posInterGen, BoltzmannCalculator bcalc) {

		this.confSpace = confSpace;
		this.posInterGen = posInterGen;
		this.bcalc = bcalc;

		staticStatic = new BigExp(1.0, 0);
		singlesPairs = new TupleMatrixGeneric<>(confSpace);
		triples = null;
		lowerBounder = false;
	}

	public ClusterZMatrix(ConfSpace confSpace, PosInterGen posInterGen, BoltzmannCalculator bcalc, boolean lowerBounder) {

		this.confSpace = confSpace;
		this.posInterGen = posInterGen;
		this.bcalc = bcalc;

		staticStatic = new BigExp(1.0, 0);
		singlesPairs = new TupleMatrixGeneric<>(confSpace);
		triples = null;
		this.lowerBounder = lowerBounder;
	}

	private BigExp z(double energy) {
		return new BigExp(bcalc.calcPrecise(energy));
	}

	public void compute(ClusterMember member, TaskExecutor tasks, boolean includeStaticStatic, Double tripleCorrectionThreshold, ConfEnergyCalculator ecalc) {

		// register with hazelcast
		member.registerService(ServiceName, this);

		if (includeStaticStatic) {
			computeStaticStatic(member, ecalc);
		}
		computeSingles(member, tasks, ecalc);
		computePairs(member, tasks, ecalc);
		if (tripleCorrectionThreshold != null) {
			computeTripleCorrections(member, tasks, ecalc, tripleCorrectionThreshold);
		}

		member.unregisterService(ServiceName);
	}

	private int calcWorkSize(TaskExecutor tasks, ConfEnergyCalculator ecalc) {
		// figure out how many singles to request at once
		// the best value here probably varies per design,
		// but a small multiple of the ecalc batch size and number of threads seems reasonable
		return ecalc.maxBatchSize()*tasks.getParallelism()*4;
	}

	private void computeStaticStatic(ClusterMember member, ConfEnergyCalculator ecalc) {

		member.log0("computing static-static ...");
		member.barrier(1, TimeUnit.MINUTES);

		// just compute the static-static energy on every node...
        if(lowerBounder){
			staticStatic = z(ecalc.calcEnergy(
					confSpace.assign(),
					posInterGen.staticStatic()
			));
		}else {
			staticStatic = z(ecalc.minimizeEnergy(
					confSpace.assign(),
					posInterGen.staticStatic()
			));
		}

		member.barrier(1, TimeUnit.MINUTES);
		member.log0("static-static finished");
	}

	private void computeSingles(ClusterMember member, TaskExecutor tasks, ConfEnergyCalculator ecalc) {

		int numSingles = singlesPairs.getNumOneBody();
		latch = new WorkLatch(numSingles);
		if (member.id() == 0) {
			progress = new Progress(numSingles);
		}
		droppedTuples = new AtomicLong(0L);
		member.log0("computing %d singles ...", numSingles);

		// make the single generator
		try (var gen = new Generator<int[]>(yielder -> {

			for (int posi=0; posi<ecalc.confSpace().numPos(); posi++) {
				for (int confi=0; confi<ecalc.confSpace().numConf(posi); confi++) {
					yielder.yield(new int[] { posi, confi });
				}
			}
		})) {

			// partiton the workload onto the cluster
			var partition = member.dynamicPartition(gen.iterator());
			member.barrier(1, TimeUnit.MINUTES);

			// do the work
			int worksSize = calcWorkSize(tasks, ecalc);
			Batch batch = new Batch(ecalc.maxBatchSize());
			while (true) {

				// get the next works
				var works = partition.nextWorks(worksSize);
				if (works.isEmpty()) {
					break;
				}

				// calculate the batches
				for (int[] work : works) {
					int posi = work[0];
					int confi = work[1];
					batch.tuples.add(new Single(posi, confi));
					if (batch.isFull()) {
						batch.submit(tasks, ecalc, member);
						batch = new Batch(ecalc.maxBatchSize());
					}
				}
				if (!batch.isEmpty()) {
					batch.submit(tasks, ecalc, member);
				}
			}
			tasks.waitForFinish();
		}

		// wait for all the other members to finish sending results
		latch.await(1, TimeUnit.MINUTES);

		// double-check that we have all the entries locally, just in case
		long missing = -droppedTuples.get();
		for (int posi=0; posi<ecalc.confSpace().numPos(); posi++) {
			for (int confi=0; confi<ecalc.confSpace().numConf(posi); confi++) {
				if (singlesPairs.getOneBody(posi, confi) == null) {
					missing += 1;
				}
			}
		}
		if (missing > 0) {
			throw new IllegalStateException("missing " + missing  + " Z matrix singles!");
		}

		// cleanup
		latch = null;
		progress = null;
		droppedTuples = null;

		member.barrier(1, TimeUnit.MINUTES);
		member.log0("singles finished");
	}

	private void computePairs(ClusterMember member, TaskExecutor tasks, ConfEnergyCalculator ecalc) {

		int numPairs = singlesPairs.getNumPairwise();
		latch = new WorkLatch(numPairs);
		if (member.id() == 0) {
			progress = new Progress(numPairs);
		}
		droppedTuples = new AtomicLong(0L);
		member.log0("computing %d pairs ...", numPairs);

		// make the pair generator
		try (var gen = new Generator<int[]>(yielder -> {

			for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
				for (int posi2=0; posi2<posi1; posi2++) {
					for (int confi1=0; confi1<ecalc.confSpace().numConf(posi1); confi1++) {
						for (int confi2=0; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
							yielder.yield(new int[] { posi1, confi1, posi2, confi2 });
						}
					}
				}
			}
		})) {

			// partiton the workload onto the cluster
			var partition = member.dynamicPartition(gen.iterator());
			member.barrier(1, TimeUnit.MINUTES);

			// do the work
			int worksSize = calcWorkSize(tasks, ecalc);
			Batch batch = new Batch(ecalc.maxBatchSize());
			while (true) {

				// get the next works
				var works = partition.nextWorks(worksSize);
				if (works.isEmpty()) {
					break;
				}

				// calculate the batches
				for (int[] work : works) {
					int posi1 = work[0];
					int confi1 = work[1];
					int posi2 = work[2];
					int confi2 = work[3];
					batch.tuples.add(new Pair(posi1, confi1, posi2, confi2));
					if (batch.isFull()) {
						batch.submit(tasks, ecalc, member);
						batch = new Batch(ecalc.maxBatchSize());
					}
				}
				if (!batch.isEmpty()) {
					batch.submit(tasks, ecalc, member);
				}
			}
			tasks.waitForFinish();
		}

		// wait for all the other members to finish sending results
		latch.await(1, TimeUnit.MINUTES);

		// double-check that we have all the entries locally, just in case
		long missing = -droppedTuples.get();
		for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				for (int confi1=0; confi1<ecalc.confSpace().numConf(posi1); confi1++) {
					for (int confi2=0; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
						if (singlesPairs.getPairwise(posi1, confi1, posi2, confi2) == null) {
							missing += 1;
						}
					}
				}
			}
		}
		if (missing > 0) {
			throw new IllegalStateException("missing " + missing + " Z matrix pairs!");
		}

		// cleanup
		latch = null;
		progress = null;
		droppedTuples = null;

		member.barrier(1, TimeUnit.MINUTES);
		member.log0("pairs finished");
	}

	private void computeTripleCorrections(ClusterMember member, TaskExecutor tasks, ConfEnergyCalculator ecalc, double energyThreshold) {

		// allocate space
		triples = new TripleMatrix<>(confSpace);

		member.log0("filtering triples ...");

		// filter the triples by the threshold
		int numTriples = 0;
		BigExp zThreshold = new BigExp(bcalc.calcPrecise(energyThreshold));
		boolean[] passed = new boolean[triples.size()];
		Arrays.fill(passed, false);
		{
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

									// triple passed the threshold!
									passed[index++] = true;
									numTriples += 1;
								}
							}
						}
					}
				}
			}
			assert (index == triples.size());
		}

		latch = new WorkLatch(numTriples);
		if (member.id() == 0) {
			progress = new Progress(numTriples);
		}
		droppedTuples = new AtomicLong(0L);
		member.log0("computing %d triples ...", numTriples);

		// make the triple generator
		try (var gen = new Generator<int[]>(yielder -> {

			int index = 0;
			for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
				for (int posi2=0; posi2<posi1; posi2++) {
					for (int posi3=0; posi3<posi2; posi3++) {
						for (int confi1=0; confi1<ecalc.confSpace().numConf(posi1); confi1++) {
							for (int confi2=0; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
								for (int confi3=0; confi3<ecalc.confSpace().numConf(posi3); confi3++) {
									if (passed[index]) {
										yielder.yield(new int[] { posi1, confi1, posi2, confi2, posi3, confi3 });
									}
									index += 1;
								}
							}
						}
					}
				}
			}
		})) {

			// partiton the workload onto the cluster
			var partition = member.dynamicPartition(gen.iterator());
			member.barrier(1, TimeUnit.MINUTES);

			// do the work
			int worksSize = calcWorkSize(tasks, ecalc);
			Batch batch = new Batch(ecalc.maxBatchSize());
			while (true) {

				// get the next works
				var works = partition.nextWorks(worksSize);
				if (works.isEmpty()) {
					break;
				}

				// calculate the batches
				for (int[] work : works) {
					int posi1 = work[0];
					int confi1 = work[1];
					int posi2 = work[2];
					int confi2 = work[3];
					int posi3 = work[4];
					int confi3 = work[5];
					batch.tuples.add(new Triple(posi1, confi1, posi2, confi2, posi3, confi3));
					if (batch.isFull()) {
						batch.submit(tasks, ecalc, member);
						batch = new Batch(ecalc.maxBatchSize());
					}
				}
				if (!batch.isEmpty()) {
					batch.submit(tasks, ecalc, member);
				}
			}
			tasks.waitForFinish();
		}

		// wait for all the other members to finish sending results
		latch.await(1, TimeUnit.MINUTES);

		// double-check that we have all the entries locally, just in case
		{
			int index = 0;
			long missing = -droppedTuples.get();
			for (int posi1=0; posi1<ecalc.confSpace().numPos(); posi1++) {
				for (int posi2=0; posi2<posi1; posi2++) {
					for (int posi3=0; posi3<posi2; posi3++) {
						for (int confi1=0; confi1<ecalc.confSpace().numConf(posi1); confi1++) {
							for (int confi2=0; confi2<ecalc.confSpace().numConf(posi2); confi2++) {
								for (int confi3=0; confi3<ecalc.confSpace().numConf(posi3); confi3++) {
									if (passed[index]) {
										if (triples.get(posi1, confi1, posi2, confi2, posi3, confi3) == null) {
											missing += 1;
										}
									}
									index += 1;
								}
							}
						}
					}
				}
			}
			if (missing > 0) {
				throw new IllegalStateException("missing " + missing + " Z matrix triple corrections!");
			}
		}

		// cleanup
		latch = null;
		progress = null;
		droppedTuples = null;

		member.barrier(1, TimeUnit.MINUTES);
		member.log0("triples finished");
	}

	private class Batch {

		final int capacity;
		final List<Tuple> tuples;

		Batch(int capacity) {
			this.capacity = capacity;
			tuples = new ArrayList<>(capacity);
		}

		int size() {
			return tuples.size();
		}

		public boolean isEmpty() {
			return size() == 0;
		}

		boolean isFull() {
			return size() == capacity;
		}

		void submit(TaskExecutor tasks, ConfEnergyCalculator ecalc, ClusterMember member) {
			assert (!isEmpty());
			tasks.submit(
				() -> {

					// compute the z values
					var jobs = tuples.stream()
						.map(tuple -> tuple.makeJob(confSpace, posInterGen))
						.collect(Collectors.toList());
					if(lowerBounder)
						ecalc.calcEnergies(jobs);
					else
						ecalc.minimizeEnergies(jobs);
					for (int i=0; i<tuples.size(); i++) {
						tuples.get(i).setZ(z(jobs.get(i).energy));
					}

					// update locally
					writeTuples(tuples);

					// update other cluster members
					member.sendToOthers(() -> new TuplesOperation(tuples));

					return 42;
				},
				answer -> {}
			);
		}
	}

	void writeTuples(List<Tuple> tuples) {
		// WARNING: this function gets hit from multiple threads!
		// We don't need to synchronize the tuple writes though.
		// Phe partition of tuples guarantees none of those writes will ever race.
		for (var tuple : tuples) {
			boolean wasWritten = tuple.write(singlesPairs, triples);
			if (!wasWritten) {
				droppedTuples.incrementAndGet();
			}
		}
		if (progress != null) {
			// we do need to prevent races on the progress bar though
			synchronized (this) {
				progress.incrementProgress(tuples.size());
			}
		}
		latch.finished(tuples.size());
	}

	public BigExp staticStatic() {
		return staticStatic;
	}

	public int numSingles() {
		return singlesPairs.getNumOneBody();
	}
	public int singleIndex(int posi, int confi) {
		return singlesPairs.getOneBodyIndex(posi, confi);
	}
	public BigExp single(int posi, int confi) {
		return singlesPairs.getOneBody(posi, confi);
	}
	public void set(int posi, int confi, BigExp val) {
		singlesPairs.setOneBody(posi, confi, val);
	}

	public int numPairs() {
		return singlesPairs.getNumPairwise();
	}
	public int pairIndex(int posi1, int posi2) {
		return singlesPairs.getPairwiseIndex(posi1, posi2);
	}
	public int pairIndex(int posi1, int confi1, int posi2, int confi2) {
		return singlesPairs.getPairwiseIndex(posi1, confi1, posi2, confi2);
	}
	public BigExp pair(int posi1, int confi1, int posi2, int confi2) {
		return singlesPairs.getPairwise(posi1, confi1, posi2, confi2);
	}
	public void set(int posi1, int confi1, int posi2, int confi2, BigExp val) {
		singlesPairs.setPairwise(posi1, confi1, posi2, confi2, val);
	}

	public boolean hasTriples() {
		return triples != null && triples.count() > 0;
	}
	public int numTriples() {
		return triples.size();
	}
	public int tripleIndex(int posi1, int posi2, int posi3) {
		return triples.index(posi1, posi2, posi3);
	}
	public int tripleIndex(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3) {
		return triples.index(posi1, confi1, posi2, confi2, posi3, confi3);
	}
	public BigExp triple(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3) {
		return triples.get(posi1, confi1, posi2, confi2, posi3, confi3);
	}
	public void set(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3, BigExp val) {
		triples.set(posi1, confi1, posi2, confi2, posi3, confi3, val);
	}

	public String singlesPairsToString(int cellWidth, int precision) {
		return singlesPairs.toString(cellWidth, f -> f.toString(precision));
	}
}
