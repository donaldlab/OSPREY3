package edu.duke.cs.osprey.coffee.seqdb;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import com.google.common.collect.Iterators;
import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.coffee.TestCoffee;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import org.junit.Test;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;


public class TestSeqDB {

	static {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();
	}

	private static final MathContext mathContext = new MathContext(16, RoundingMode.HALF_UP);
	private static final BigDecimalBounds emptySum = BigDecimalBounds.makeZero();
	private static final BigDecimalBounds unknownBound = new BigDecimalBounds(BigDecimal.ZERO, MathTools.BigPositiveInfinity);

	private static void withSeqDB(MultiStateConfSpace confSpace, Consumer<SeqDB> block) {
		withSeqDBs(confSpace, 1, block);
	}

	private static void withSeqDBs(MultiStateConfSpace confSpace, int numMembers, Consumer<SeqDB> block) {
		var exceptions = ClusterMember.launchPseudoCluster(numMembers, cluster -> {
			try (var member = new ClusterMember(cluster)) {

				// make the sequence database
				try (var seqdb = new SeqDB(confSpace, mathContext, null, member)) {

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
			assertThat(seqdb.boundsUnsequenced(target), is(unknownBound));

			assertThat(Iterators.size(seqdb.getSums().iterator()), is(0));
			assertThat(Iterators.size(seqdb.boundsSequenced().iterator()), is(0));

			// check partial sequence
			var seq = confSpace.seqSpace.makeUnassignedSequence();
			assertThat(seqdb.getSums(seq), is(nullValue()));
			assertThat(seqdb.boundsUnexplored(seq).get(complex), is(unknownBound));
			assertThat(seqdb.boundsUnexplored(seq).get(design), is(unknownBound));

			// check full sequence
			seq = confSpace.seqSpace.makeUnassignedSequence()
				.set(confSpace.seqSpace.positions.get(0), confSpace.seqSpace.positions.get(0).resTypes.get(0));
			assertThat(seqdb.bounds(seq).get(complex), is(unknownBound));
			assertThat(seqdb.bounds(seq).get(design), is(unknownBound));
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
			assertThat(seqdb.getSums(seq).get(complex), is(new BigDecimalBounds(0.0, 3.0)));
			assertThat(seqdb.getSums(seq).get(design), is(new BigDecimalBounds(0.0, 5.0)));
			assertThat(seqdb.boundsUnsequenced(target), is(new BigDecimalBounds(0.0, 7.0)));
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
				seqdb.member.waitForOperation(ops + 1, 2, TimeUnit.SECONDS);
			}
			seqdb.member.barrier(2, TimeUnit.SECONDS);

			if (seqdb.member.isDirector()) {
				// should have been added
				assertThat(seqdb.getSums(seq).get(complex), is(new BigDecimalBounds(0.0, 3.0)));
				assertThat(seqdb.getSums(seq).get(design), is(new BigDecimalBounds(0.0, 5.0)));
				assertThat(seqdb.boundsUnsequenced(target), is(new BigDecimalBounds(0.0, 7.0)));
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
			assertThat(seqdb.getSums(seq).get(complex), is(new BigDecimalBounds(0.0, 12.0)));
			assertThat(seqdb.getSums(seq).get(design), is(new BigDecimalBounds(0.0, 13.0)));
			assertThat(seqdb.boundsUnsequenced(target), is(new BigDecimalBounds(0.0, 14.0)));
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
			assertThat(seqdb.getSums(seq).get(complex), is(new BigDecimalBounds(0.0, 2.0)));
			assertThat(seqdb.getSums(seq).get(design), is(new BigDecimalBounds(0.0, 3.0)));
			assertThat(seqdb.boundsUnsequenced(target), is(new BigDecimalBounds(0.0, 4.0)));
		});
	}
}
