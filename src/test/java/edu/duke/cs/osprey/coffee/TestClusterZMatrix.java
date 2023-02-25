package edu.duke.cs.osprey.coffee;

import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.*;
import static org.junit.jupiter.api.Assertions.fail;

import edu.duke.cs.osprey.coffee.zmat.ClusterZMatrix;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.junit.jupiter.api.Test;

import java.util.function.Consumer;




public class TestClusterZMatrix {

	static {
		Cluster.fixHazelcastLogging();
	}

	@Test
	public void pos3_desmet_local() {
		compareZMat(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"complex",
			PosInterDist.DesmetEtAl1992,
			null,
			0,
			4
		);
	}

	@Test
	public void pos3_desmet_1member() {
		compareZMat(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"complex",
			PosInterDist.DesmetEtAl1992,
			null,
			1,
			4
		);
	}

	@Test
	public void pos3_desmet_2members() {
		compareZMat(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"complex",
			PosInterDist.DesmetEtAl1992,
			null,
			2,
			2
		);
	}

	@Test
	public void pos3_desmet_4members() {
		compareZMat(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"complex",
			PosInterDist.DesmetEtAl1992,
			null,
			4,
			1
		);
	}

	@Test
	public void pos6_desmet_local() {
		compareZMat(
			TestCoffee.affinity_6ov7_2mut4flex(),
			"complex",
			PosInterDist.DesmetEtAl1992,
			10.0,
			0,
			4
		);
	}

	@Test
	public void pos6_desmet_4members() {
		compareZMat(
			TestCoffee.affinity_6ov7_2mut4flex(),
			"complex",
			PosInterDist.DesmetEtAl1992,
			10.0,
			4,
			1
		);
	}

	/** test that the clustered version computes the same matrix as the local version */
	private void compareZMat(MultiStateConfSpace confSpace, String stateName, PosInterDist posInterDist, Double triplesThreshold, int numMembers, int numThreads) {

		// compute a reference zmat locally
		var zmatExp = calcZMat(null, confSpace, stateName, posInterDist, triplesThreshold, Math.max(numMembers, 1)*numThreads);

		if (numMembers <= 0) {
			return;
		}

		withPseudoCluster(numMembers, cluster -> {

			var zmatObs = calcZMat(cluster, confSpace, stateName, posInterDist, triplesThreshold, numThreads);

			// compare the zmats
			var cs = confSpace.getState(stateName).confSpace;
			for (int posi1=0; posi1<cs.numPos(); posi1++) {
				for (int posi2=0; posi2<posi1; posi2++) {
					for (int confi1=0; confi1<cs.numConf(posi1); confi1++) {
						assertThat(zmatObs.single(posi1, confi1), is(zmatExp.single(posi1, confi1)));
						for (int confi2=0; confi2<cs.numConf(posi2); confi2++) {
							assertThat(zmatObs.pair(posi1, confi1, posi2, confi2), is(zmatExp.pair(posi1, confi1, posi2, confi2)));
						}
					}
				}
			}
		});
	}

	private static void withPseudoCluster(int numMembers, Consumer<Cluster> block) {
		var exceptions = ClusterMember.launchPseudoCluster(numMembers, block);
		if (!exceptions.isEmpty()) {
			fail("Cluster threads encountered exceptions");
		}
	}

	private ClusterZMatrix calcZMat(Cluster cluster, MultiStateConfSpace confSpace, String stateName, PosInterDist posInterDist, Double triplesThreshold, int numThreads) {

		var coffee = new Coffee.Builder(confSpace)
			.setCluster(cluster)
			.setParallelism(Parallelism.makeCpu(numThreads))
			.configEachState(config -> {
				config.posInterGen = new PosInterGen(posInterDist, null);
			})
			.setStaticStatic(false)
			.setTripleCorrectionThreshold(triplesThreshold)
			.build();

		var state = confSpace.getState(stateName);

		// compute the zmat
		return coffee.calcZMat(state.index);
	}
}
