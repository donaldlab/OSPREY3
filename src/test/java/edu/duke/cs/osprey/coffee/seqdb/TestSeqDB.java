package edu.duke.cs.osprey.coffee.seqdb;

import static edu.duke.cs.osprey.TestBase.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import com.google.common.collect.Iterators;
import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.coffee.TestCoffee;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import org.junit.Test;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.math.BigDecimal;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestSeqDB {

	static {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();
	}

	private static final BigDecimalBounds unknownBound = new BigDecimalBounds(BigDecimal.ZERO, MathTools.BigPositiveInfinity);

	private static void withSeqDB(MultiStateConfSpace confSpace, Consumer<SeqDB> block) {
		withSeqDBs(confSpace, 1, block);
	}

	private static void withSeqDBs(MultiStateConfSpace confSpace, int numMembers, Consumer<SeqDB> block) {
		var exceptions = ClusterMember.launchPseudoCluster(numMembers, cluster -> {
			try (var member = new ClusterMember(cluster)) {

				// make the sequence database
				try (var seqdb = new SeqDB.Builder(confSpace, member)
					.setNumBestConfs(10)
					.build()
				) {

					// wait for all the database instances to be ready
					member.barrier(1, TimeUnit.MINUTES);

					block.accept(seqdb);
				}
			}
		});
		if (!exceptions.isEmpty()) {
			fail("Cluster threads encountered exceptions");
		}
	}

	@Test
	public void empty() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		withSeqDB(confSpace, seqdb -> {

			// there shouldn't be any info by default

			// check unsequenced states
			assertThat(seqdb.getSum(target), is(nullValue()));
			assertThat(seqdb.getUnsequenced(target).zSumBounds, is(unknownBound));
			assertThat(seqdb.getUnsequenced(target).zSumDropped, is(BigDecimal.ZERO));

			assertThat(Iterators.size(seqdb.getSums().iterator()), is(0));
			assertThat(Iterators.size(seqdb.getSequenced().iterator()), is(0));

			// check partial sequence
			var seq = confSpace.seqSpace.makeUnassignedSequence();
			assertThat(seqdb.getSums(seq), is(nullValue()));
			assertThat(seqdb.getUnexplored(seq).get(complex).zSumBounds, is(unknownBound));
			assertThat(seqdb.getUnexplored(seq).get(design).zSumBounds, is(unknownBound));
			assertThat(seqdb.getUnexplored(seq).get(complex).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getUnexplored(seq).get(design).zSumDropped, is(BigDecimal.ZERO));

			// check full sequence
			seq = confSpace.seqSpace.makeUnassignedSequence()
				.set(confSpace.seqSpace.positions.get(0), confSpace.seqSpace.positions.get(0).resTypes.get(0));
			assertThat(seqdb.get(seq).get(complex).zSumBounds, is(unknownBound));
			assertThat(seqdb.get(seq).get(design).zSumBounds, is(unknownBound));
			assertThat(seqdb.get(seq).get(complex).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.get(seq).get(design).zSumDropped, is(BigDecimal.ZERO));
		});
	}

	@Test
	public void saveEmpty() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State target = confSpace.getState("target");

		withSeqDB(confSpace, seqdb -> {

			var batch = seqdb.batch();
			batch.save();

			// there still shouldn't be any info
			assertThat(seqdb.getSum(target), is(nullValue()));
			assertThat(Iterators.size(seqdb.getSums().iterator()), is(0));
		});
	}

	@Test
	public void addLocal() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		withSeqDB(confSpace, seqdb -> {

			// add one value for each state
			var batch = seqdb.batch();
			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();
			batch.addZSumUpper(complex, seq, new BigExp(3.0));
			batch.addZSumUpper(design, seq, new BigExp(5.0));
			batch.addZSumUpper(target, null, new BigExp(7.0));
			batch.save();

			// all values should have been added
			assertThat(seqdb.getSums(seq).get(complex).zSumBounds, is(new BigDecimalBounds(0.0, 3.0)));
			assertThat(seqdb.getSums(seq).get(design).zSumBounds, is(new BigDecimalBounds(0.0, 5.0)));
			assertThat(seqdb.getSum(target).zSumBounds, is(new BigDecimalBounds(0.0, 7.0)));

			// nothing should have been dropped
			assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSum(target).zSumDropped, is(BigDecimal.ZERO));
		});
	}

	@Test
	public void addRemote() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		withSeqDBs(confSpace, 2, seqdb -> {

			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();
			
			long ops = seqdb.member.finishedOperations();

			if (!seqdb.member.isDirector()) {
				// add one value for each state
				var batch = seqdb.batch();
				batch.addZSumUpper(complex, seq, new BigExp(3.0));
				batch.addZSumUpper(design, seq, new BigExp(5.0));
				batch.addZSumUpper(target, null, new BigExp(7.0));
				batch.save();
			}

			// wait for the save to finish
			if (seqdb.member.isDirector()) {
				seqdb.member.waitForOperation(ops + 1, 5, TimeUnit.SECONDS);
			}
			seqdb.member.barrier(10, TimeUnit.SECONDS);

			if (seqdb.member.isDirector()) {

				// everything should have been added
				assertThat(seqdb.getSums(seq).get(complex).zSumBounds, is(new BigDecimalBounds(0.0, 3.0)));
				assertThat(seqdb.getSums(seq).get(design).zSumBounds, is(new BigDecimalBounds(0.0, 5.0)));
				assertThat(seqdb.getSum(target).zSumBounds, is(new BigDecimalBounds(0.0, 7.0)));

				// nothing should have been dropped
				assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(BigDecimal.ZERO));
				assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(BigDecimal.ZERO));
				assertThat(seqdb.getSum(target).zSumDropped, is(BigDecimal.ZERO));
			}

			seqdb.member.barrier(2, TimeUnit.SECONDS);
		});
	}

	@Test
	public void addConfsLocal() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		assertThat(complex.confSpace.numPos(), is(3));
		assertThat(design.confSpace.numPos(), is(1));
		assertThat(target.confSpace.numPos(), is(2));

		withSeqDB(confSpace, seqdb -> {

			var confComplex = new ConfSearch.EnergiedConf(new int[] { 1, 2, 3 }, -4.0, -3.0);
			var confDesign = new ConfSearch.EnergiedConf(new int[] { 1 }, -6.0, -5.0);
			var confTarget = new ConfSearch.EnergiedConf(new int[] { 1, 2 }, -8.0, -7.0);

			// add one value for each state
			var batch = seqdb.batch();
			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();
			batch.addZConf(complex, seq, new BigDecimal("3.0"), new BigExp(0.0), confComplex);
			batch.addZConf(design, seq, new BigDecimal("5.0"), new BigExp(0.0), confDesign);
			batch.addZConf(target, null, new BigDecimal("7.0"), new BigExp(0.0), confTarget);
			batch.save();

			// all values should have been added
			assertThat(seqdb.getSums(seq).get(complex).zSumBounds, is(new BigDecimalBounds(3.0, 3.0)));
			assertThat(seqdb.getSums(seq).get(design).zSumBounds, is(new BigDecimalBounds(5.0, 5.0)));
			assertThat(seqdb.getSum(target).zSumBounds, is(new BigDecimalBounds(7.0, 7.0)));

			// nothing should have been dropped
			assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSum(target).zSumDropped, is(BigDecimal.ZERO));

			// check best confs
			assertThat(seqdb.getBestConfs(complex, seq), contains(confComplex));
			assertThat(seqdb.getBestConfs(design, seq), contains(confDesign));
			assertThat(seqdb.getBestConfs(target), contains(confTarget));
		});
	}

	@Test
	public void addAddConfsLocal() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		assertThat(complex.confSpace.numPos(), is(3));
		assertThat(design.confSpace.numPos(), is(1));
		assertThat(target.confSpace.numPos(), is(2));

		withSeqDB(confSpace, seqdb -> {

			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();
			ConfSearch.EnergiedConf[] confComplex = {
				new ConfSearch.EnergiedConf(new int[] { 1, 2, 3 }, -4.0, -3.0),
				new ConfSearch.EnergiedConf(new int[] { 1, 2, 0 }, -6.0, -5.0)
			};
			ConfSearch.EnergiedConf[] confDesign = {
				new ConfSearch.EnergiedConf(new int[] { 1 }, -6.0, -5.0),
				new ConfSearch.EnergiedConf(new int[] { 0 }, -8.0, -7.0)
			};
			ConfSearch.EnergiedConf[] confTarget = {
				new ConfSearch.EnergiedConf(new int[] { 1, 2 }, -8.0, -7.0),
				new ConfSearch.EnergiedConf(new int[] { 1, 0 }, -10.0, -9.0)
			};

			// add one value for each state
			var batch = seqdb.batch();
			batch.addZConf(complex, seq, new BigDecimal("3.0"), new BigExp(0.0), confComplex[0]);
			batch.addZConf(design, seq, new BigDecimal("5.0"), new BigExp(0.0), confDesign[0]);
			batch.addZConf(target, null, new BigDecimal("7.0"), new BigExp(0.0), confTarget[0]);
			batch.save();

			// add one value for each state
			batch = seqdb.batch();
			batch.addZConf(complex, seq, new BigDecimal("1.0"), new BigExp(0.0), confComplex[1]);
			batch.addZConf(design, seq, new BigDecimal("2.0"), new BigExp(0.0), confDesign[1]);
			batch.addZConf(target, null, new BigDecimal("3.0"), new BigExp(0.0), confTarget[1]);
			batch.save();

			// all values should have been added
			assertThat(seqdb.getSums(seq).get(complex).zSumBounds, is(new BigDecimalBounds(4.0, 4.0)));
			assertThat(seqdb.getSums(seq).get(design).zSumBounds, is(new BigDecimalBounds(7.0, 7.0)));
			assertThat(seqdb.getSum(target).zSumBounds, is(new BigDecimalBounds(10.0, 10.0)));

			// nothing should have been dropped
			assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSum(target).zSumDropped, is(BigDecimal.ZERO));

			// check best confs
			assertThat(seqdb.getBestConfs(complex, seq), contains(confComplex[1], confComplex[0]));
			assertThat(seqdb.getBestConfs(design, seq), contains(confDesign[1], confDesign[0]));
			assertThat(seqdb.getBestConfs(target), contains(confTarget[1], confTarget[0]));
		});
	}

	@Test
	public void addConfsRemote() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		assertThat(complex.confSpace.numPos(), is(3));
		assertThat(design.confSpace.numPos(), is(1));
		assertThat(target.confSpace.numPos(), is(2));

		withSeqDBs(confSpace, 2, seqdb -> {

			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();
			var confComplex = new ConfSearch.EnergiedConf(new int[] { 1, 2, 3 }, -4.0, -3.0);
			var confDesign = new ConfSearch.EnergiedConf(new int[] { 1 }, -6.0, -5.0);
			var confTarget = new ConfSearch.EnergiedConf(new int[] { 1, 2 }, -8.0, -7.0);

			long ops = seqdb.member.finishedOperations();

			if (!seqdb.member.isDirector()) {
				// add one value for each state
				var batch = seqdb.batch();
				batch.addZConf(complex, seq, new BigDecimal("3.0"), new BigExp(0.0), confComplex);
				batch.addZConf(design, seq, new BigDecimal("5.0"), new BigExp(0.0), confDesign);
				batch.addZConf(target, null, new BigDecimal("7.0"), new BigExp(0.0), confTarget);
				batch.save();
			}

			// wait for the save to finish
			if (seqdb.member.isDirector()) {
				seqdb.member.waitForOperation(ops + 1, 5, TimeUnit.SECONDS);
			}
			seqdb.member.barrier(10, TimeUnit.SECONDS);

			if (seqdb.member.isDirector()) {

				// all values should have been added
				assertThat(seqdb.getSums(seq).get(complex).zSumBounds, is(new BigDecimalBounds(3.0, 3.0)));
				assertThat(seqdb.getSums(seq).get(design).zSumBounds, is(new BigDecimalBounds(5.0, 5.0)));
				assertThat(seqdb.getSum(target).zSumBounds, is(new BigDecimalBounds(7.0, 7.0)));

				// nothing should have been dropped
				assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(BigDecimal.ZERO));
				assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(BigDecimal.ZERO));
				assertThat(seqdb.getSum(target).zSumDropped, is(BigDecimal.ZERO));

				// check best confs
				assertThat(seqdb.getBestConfs(complex, seq), contains(confComplex));
				assertThat(seqdb.getBestConfs(design, seq), contains(confDesign));
				assertThat(seqdb.getBestConfs(target), contains(confTarget));
			}

			seqdb.member.barrier(2, TimeUnit.SECONDS);
		});
	}

	@Test
	public void addAddLocal() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		withSeqDB(confSpace, seqdb -> {

			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();

			// add one value for each state
			var batch = seqdb.batch();
			batch.addZSumUpper(complex, seq, new BigExp(3.0));
			batch.addZSumUpper(design, seq, new BigExp(5.0));
			batch.addZSumUpper(target, null, new BigExp(7.0));
			batch.save();

			// do it again
			batch = seqdb.batch();
			batch.addZSumUpper(complex, seq, new BigExp(9.0));
			batch.addZSumUpper(design, seq, new BigExp(8.0));
			batch.addZSumUpper(target, null, new BigExp(7.0));
			batch.save();

			// all values should have been added together
			assertThat(seqdb.getSums(seq).get(complex).zSumBounds, is(new BigDecimalBounds(0.0, 12.0)));
			assertThat(seqdb.getSums(seq).get(design).zSumBounds, is(new BigDecimalBounds(0.0, 13.0)));
			assertThat(seqdb.getSum(target).zSumBounds, is(new BigDecimalBounds(0.0, 14.0)));

			// nothing should have been dropped
			assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSum(target).zSumDropped, is(BigDecimal.ZERO));
		});
	}

	@Test
	public void addSubLocal() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		withSeqDB(confSpace, seqdb -> {

			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();

			// add one value for each state
			var batch = seqdb.batch();
			batch.addZSumUpper(complex, seq, new BigExp(3.0));
			batch.addZSumUpper(design, seq, new BigExp(5.0));
			batch.addZSumUpper(target, null, new BigExp(7.0));
			batch.save();

			// do it again
			batch = seqdb.batch();
			batch.subZSumUpper(complex, seq, new BigExp(1.0));
			batch.subZSumUpper(design, seq, new BigExp(2.0));
			batch.subZSumUpper(target, null, new BigExp(3.0));
			batch.save();

			// all values should have been added together
			assertThat(seqdb.getSums(seq).get(complex).zSumBounds, is(new BigDecimalBounds(0.0, 2.0)));
			assertThat(seqdb.getSums(seq).get(design).zSumBounds, is(new BigDecimalBounds(0.0, 3.0)));
			assertThat(seqdb.getSum(target).zSumBounds, is(new BigDecimalBounds(0.0, 4.0)));

			// nothing should have been dropped
			assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(BigDecimal.ZERO));
			assertThat(seqdb.getSum(target).zSumDropped, is(BigDecimal.ZERO));
		});
	}

	@Test
	public void dropLocal() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		withSeqDB(confSpace, seqdb -> {

			// drop one for each state
			var batch = seqdb.batch();
			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();
			batch.drop(complex, seq, new BigExp(3.0));
			batch.drop(design, seq, new BigExp(5.0));
			batch.drop(target, null, new BigExp(7.0));
			batch.save();

			// check drops
			assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(MathTools.biggen(3.0)));
			assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(MathTools.biggen(5.0)));
			assertThat(seqdb.getSum(target).zSumDropped, is(MathTools.biggen(7.0)));
		});
	}

	@Test
	public void dropRemote() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		withSeqDBs(confSpace, 2, seqdb -> {

			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();

			long ops = seqdb.member.finishedOperations();

			if (!seqdb.member.isDirector()) {
				// drop one for each state
				var batch = seqdb.batch();
				batch.drop(complex, seq, new BigExp(3.0));
				batch.drop(design, seq, new BigExp(5.0));
				batch.drop(target, null, new BigExp(7.0));
				batch.save();
			}

			// wait for the save to finish
			if (seqdb.member.isDirector()) {
				seqdb.member.waitForOperation(ops + 1, 5, TimeUnit.SECONDS);
			}
			seqdb.member.barrier(10, TimeUnit.SECONDS);

			if (seqdb.member.isDirector()) {

				// check drops
				assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(MathTools.biggen(3.0)));
				assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(MathTools.biggen(5.0)));
				assertThat(seqdb.getSum(target).zSumDropped, is(MathTools.biggen(7.0)));
			}

			seqdb.member.barrier(2, TimeUnit.SECONDS);
		});
	}

	@Test
	public void dropDropLocal() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
		MultiStateConfSpace.State complex = confSpace.getState("complex");
		MultiStateConfSpace.State design = confSpace.getState("design");
		MultiStateConfSpace.State target = confSpace.getState("target");

		withSeqDB(confSpace, seqdb -> {

			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();

			// drop one for each state
			var batch = seqdb.batch();
			batch.drop(complex, seq, new BigExp(3.0));
			batch.drop(design, seq, new BigExp(5.0));
			batch.drop(target, null, new BigExp(7.0));
			batch.save();

			// do it again
			batch = seqdb.batch();
			batch.drop(complex, seq, new BigExp(9.0));
			batch.drop(design, seq, new BigExp(8.0));
			batch.drop(target, null, new BigExp(7.0));
			batch.save();

			// check drops
			assertThat(seqdb.getSums(seq).get(complex).zSumDropped, is(MathTools.biggen(12.0)));
			assertThat(seqdb.getSums(seq).get(design).zSumDropped, is(MathTools.biggen(13.0)));
			assertThat(seqdb.getSum(target).zSumDropped, is(MathTools.biggen(14.0)));
		});
	}

	private void addLotsLocal(int numThreads) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();

		withSeqDB(confSpace, seqdb -> {

			var state = confSpace.states.get(0);
			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();

			// generate a bunch of random values
			var rand = new Random(12345);
			List<BigExp> allValues = IntStream.range(0, 1_000_000)
				.mapToObj(i -> new BigExp(rand.nextDouble(), rand.nextInt(100)))
				.collect(Collectors.toList());

			// add them up
			var bm = seqdb.bigMath()
				.set(0.0);
			for (var val : allValues) {
				bm.add(val);
			}
			var sum = bm.get();

			// add the values from a bunch of threads
			List<Thread> threads = IntStream.range(0, numThreads)
				.mapToObj(t -> new Thread(() -> {

					Batch batch = null;
					int batchSize = 0;

					// add a slice of the values to the db, in small batches
					int size = MathTools.divUp(allValues.size(), numThreads);
					int start = t*size;
					int stop = Math.min(start + size, allValues.size());
					for (var val : allValues.subList(start, stop)) {

						if (batch == null) {
							batch = seqdb.batch();
						}

						batch.addZSumUpper(state, seq, val);

						batchSize++;
						if (batchSize == 10) {
							batch.save();
							batch = null;
							batchSize = 0;
						}
					}
					if (batch != null) {
						batch.save();
					}
				}))
				.collect(Collectors.toList());
			threads.forEach(t -> t.start());
			threads.forEach(t -> {
				try {
					t.join();
				} catch (InterruptedException ex) {
					throw new RuntimeException(ex);
				}
			});

			// check the sum
			assertThat(seqdb.getSums(seq).get(state).zSumBounds.upper, isRelatively(sum, 1e-6));
		});
	}
	@Test public void addLotsLocal_1() { addLotsLocal(1); }
	@Test public void addLotsLocal_2() { addLotsLocal(2); }
	@Test public void addLotsLocal_4() { addLotsLocal(4); }
}
