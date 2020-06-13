package edu.duke.cs.osprey.coffee;

import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.*;
import static edu.duke.cs.osprey.tools.Log.log;

import edu.duke.cs.osprey.coffee.directors.SequenceDirector;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Progress;
import org.junit.Test;

import java.util.Arrays;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestCoffee {

	private static final Parallelism oneCpu = Parallelism.makeCpu(1);
	private static final Parallelism allCpus = Parallelism.makeCpu(Parallelism.getMaxNumCPUs());

	private static final double freeEnergyEpsilon = 1e-6;

	/**
	 * A 7-residue design with mutable residues on both sides
	 */
	public static MultiStateConfSpace affinity_2RL0_7mut() {
		var confSpaces = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();
		return new MultiStateConfSpace
			.Builder("complex", confSpaces.complex)
			.addMutableState("A", confSpaces.chainA)
			.addMutableState("B", confSpaces.chainB)
			.build();
	}

	/**
	 * A tiny peptide inhibitor affinity design
	 */
	public static MultiStateConfSpace affinity_6ov7_1mut2flex() {
		var confSpaces = new TestConfSpace.AffinityCompiled(
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.tiny.complex.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.tiny.design.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.tiny.target.ccsx"))
		);
		return new MultiStateConfSpace
			.Builder("complex", confSpaces.complex)
			.addMutableState("design", confSpaces.chainA)
			.addUnmutableState("target", confSpaces.chainB)
			.build();
	}

	private static void bruteForce_affinity_6ov7_1mut2flex() {
		var coffee = makeCoffee(TestCoffee.affinity_6ov7_1mut2flex(), null, allCpus);
		for (var seq : coffee.confSpace.seqSpace.getSequences()) {
			bruteForceFreeEnergies(coffee, seq);
		}

		//sequence [6 GLN=gln]
		//	state complex = -1377.127950
		//	state design = -144.199934
		//	state target = -1187.667391
		//sequence [6 GLN=ALA]
		//	state complex = -1363.561940
		//	state design = -132.431356
		//	state target = -1187.667391
		//sequence [6 GLN=ASN]
		//	state complex = -1375.773406
		//	state design = -143.920583
		//	state target = -1187.667391
	}

	public static void main(String[] args) {
		bruteForce_affinity_6ov7_1mut2flex();
	}

	private static void withPseudoCluster(int numMembers, Consumer<Cluster> block) {
		var exceptions = ClusterMember.launchPseudoCluster(numMembers, block);
		if (!exceptions.isEmpty()) {
			fail("Cluster threads encountered exceptions");
		}
	}

	private static Coffee makeCoffee(MultiStateConfSpace confSpace, Cluster cluster, Parallelism parallelism) {

		// looser bounds makes the zmat calculation much faster, and the overall computation faster for small conf spaces
		var posInterDist = PosInterDist.DesmetEtAl1992;
		//var posInterDist = PosInterDist.TighterBounds;

		return new Coffee.Builder(confSpace)
			.setCluster(cluster)
			.setParallelism(parallelism)
			.configEachState(config -> {
				config.ecalc = new CPUConfEnergyCalculator((ConfSpace)config.state.confSpace);
				config.posInterGen = new PosInterGen(posInterDist, null);
			})
			.build();
	}

	private static void bruteForceFreeEnergies(Coffee coffee, Sequence seq) {

		log("sequence [%s]", seq);

		var gcalc = new FreeEnergyCalculator();

		try (var tasks = new ThreadPoolTaskExecutor()) {
			tasks.start(coffee.parallelism.numThreads);

			for (var state : coffee.confSpace.states) {
				var stateInfo = coffee.infos[state.index];

				// collect all the conformations in one big iterator
				var rcs = seq.makeRCs(state.confSpace);
				var allConfs = MathTools.cartesianProduct(
					IntStream.range(0, state.confSpace.numPos())
						.mapToObj(posi -> Arrays.stream(rcs.get(posi))
							.boxed()
							.collect(Collectors.toList())
						)
						.collect(Collectors.toList())
				);

				// calc all the boltzmann-weighted energies and sum them together
				var bigMath = new BigMath(coffee.mathContext)
					.set(0);
				var progress = new Progress(rcs.getNumConformations().longValueExact());
				for (var confList : allConfs) {
					tasks.submit(
						() -> {
							var conf = confList.stream()
								.mapToInt(i -> i)
								.toArray();
							return stateInfo.zPath(conf);
						},
						z -> {
							bigMath.add(z);
							progress.incrementProgress();
						}
					);
				}
				tasks.waitForFinish();

				// finally, calculate the free energy
				double g = gcalc.calc(bigMath.get());
				log("\tstate %s = %f", state.name, g);
			}
		}
	}

	public void seqFreeEnergy(MultiStateConfSpace confSpace, Function<SeqSpace,Sequence> seqFunc, double precision, double[] freeEnergies, int numMembers, Parallelism parallelism) {
		withPseudoCluster(numMembers, cluster -> {

			// get the sequence
			Coffee coffee = makeCoffee(confSpace, cluster, parallelism);
			Sequence seq = seqFunc.apply(coffee.confSpace.seqSpace);

			// run COFFEE
			var director = new SequenceDirector(coffee.confSpace, seq, precision);
			coffee.run(director);

			// check the free energies
			if (cluster.nodeId == 0) {
				for (var state : confSpace.states) {
					assertThat(state.name, director.getFreeEnergy(state), isAbsoluteBound(freeEnergies[state.index], freeEnergyEpsilon));
				}
			}
		});
	}


	// the basic test
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x1() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			0.1, new double[] { -1377.127950, -144.199934, -1187.667391 },
			1, oneCpu
		);
	}

	// vary precision
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_001_1x1() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			0.01, new double[] { -1377.127950, -144.199934, -1187.667391 },
			1, oneCpu
		);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_0001_1x1() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			0.001, new double[] { -1377.127950, -144.199934, -1187.667391 },
			1, oneCpu
		);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_0_1x1() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			0.0001, new double[] { -1377.127950, -144.199934, -1187.667391 },
			1, oneCpu
		);
	}

	// vary sequence
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_ala_01_1x1() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ALA"),
			0.1, new double[] { -1363.561940, -132.431356, -1187.667391 },
			1, oneCpu
		);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_asn_01_1x1() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ASN"),
			0.1, new double[] { -1375.773406, -143.920583, -1187.667391 },
			1, oneCpu
		);
	}

	// vary cluster members/parallelism
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_2x1() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			0.1, new double[] { -1377.127950, -144.199934, -1187.667391 },
			2, oneCpu
		);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_2x2() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			0.1, new double[] { -1377.127950, -144.199934, -1187.667391 },
			2, Parallelism.makeCpu(2)
		);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x2() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			0.1, new double[] { -1377.127950, -144.199934, -1187.667391 },
			1, Parallelism.makeCpu(2)
		);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x4() {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			0.1, new double[] { -1377.127950, -144.199934, -1187.667391 },
			1, Parallelism.makeCpu(4)
		);
	}
}
