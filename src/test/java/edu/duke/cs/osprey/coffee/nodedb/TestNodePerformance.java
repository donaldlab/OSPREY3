package edu.duke.cs.osprey.coffee.nodedb;

import static edu.duke.cs.osprey.TestBase.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.coffee.TestCoffee;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import org.junit.Test;


public class TestNodePerformance {

	private static final MultiStateConfSpace confSpace = TestCoffee.affinity_6ov7_1mut2flex();
	private static final int statei = 0;
	private static final int[] conf = { -1, -1, -1 };
	private static final BigExp zero = new BigExp(0.0, 0);
	private static final BigExp zSumUpper = new BigExp(1.0, 10);

	@Test
	public void beforeUpdates() {

		var perf = new NodePerformance(confSpace);
		var score = perf.score(statei, conf, zSumUpper);

		// the score should be something reasonable
		assertThat(score.greaterThan(zero), is(true));
		assertThat(score.isFinite(), is(true));
	}

	@Test
	public void afterOneUpdateSame() {

		var perf = new NodePerformance(confSpace);
		var initialScore = perf.score(statei, conf, zSumUpper);
		perf.update(statei, conf, NodePerformance.InitialNs);
		var score = perf.score(statei, conf, zSumUpper);

		// the score should be exactly the same as the initial score
		assertThat(score.greaterThan(zero), is(true));
		assertThat(score.isFinite(), is(true));
		assertThat(score, isRelatively(initialScore, 1e-6));
	}

	@Test
	public void afterManyUpdatesSame() {

		var perf = new NodePerformance(confSpace);
		var initialScore = perf.score(statei, conf, zSumUpper);
		for (int i=0; i<NodePerformance.HistorySize; i++) {
			perf.update(statei, conf, NodePerformance.InitialNs);
		}
		var score = perf.score(statei, conf, zSumUpper);

		// the score should be exactly the same as the initial score
		assertThat(score.greaterThan(zero), is(true));
		assertThat(score.isFinite(), is(true));
		assertThat(score, isRelatively(initialScore, 1e-6));
	}

	@Test
	public void afterOneUpdateSlower() {

		var perf = new NodePerformance(confSpace);
		var initialScore = perf.score(statei, conf, zSumUpper);
		perf.update(statei, conf, NodePerformance.InitialNs*10);
		var score = perf.score(statei, conf, zSumUpper);

		// the score should still be close to the initial score, but less
		assertThat(score.greaterThan(zero), is(true));
		assertThat(score.isFinite(), is(true));
		assertThat(score.lessThan(initialScore), is(true));
		assertThat(score, isRelatively(initialScore, 1e-1));
	}

	@Test
	public void afterManyUpdatesSlower() {

		var perf = new NodePerformance(confSpace);
		var initialScore = perf.score(statei, conf, zSumUpper);
		for (int i=0; i<NodePerformance.HistorySize; i++) {
			perf.update(statei, conf, NodePerformance.InitialNs*10);
		}
		var score = perf.score(statei, conf, zSumUpper);

		// the score should be less than the initial score
		assertThat(score.greaterThan(zero), is(true));
		assertThat(score.isFinite(), is(true));
		assertThat(score.lessThan(initialScore), is(true));
	}

	@Test
	public void afterOneUpdateFaster() {

		var perf = new NodePerformance(confSpace);
		var initialScore = perf.score(statei, conf, zSumUpper);
		perf.update(statei, conf, NodePerformance.InitialNs/10);
		var score = perf.score(statei, conf, zSumUpper);

		// the score should still be close to the initial score, but greater
		assertThat(score.greaterThan(zero), is(true));
		assertThat(score.isFinite(), is(true));
		assertThat(score.greaterThan(initialScore), is(true));
		assertThat(score, isRelatively(initialScore, 1e-1));
	}

	@Test
	public void afterManyUpdatesFaster() {

		var perf = new NodePerformance(confSpace);
		var initialScore = perf.score(statei, conf, zSumUpper);
		for (int i=0; i<NodePerformance.HistorySize; i++) {
			perf.update(statei, conf, NodePerformance.InitialNs/10);
		}
		var score = perf.score(statei, conf, zSumUpper);

		// the score should be greater than the initial score
		assertThat(score.greaterThan(zero), is(true));
		assertThat(score.isFinite(), is(true));
		assertThat(score.greaterThan(initialScore), is(true));
	}
}
