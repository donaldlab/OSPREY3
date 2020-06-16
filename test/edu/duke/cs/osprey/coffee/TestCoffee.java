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
import edu.duke.cs.osprey.energy.compiled.CudaConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Progress;
import org.junit.Test;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestCoffee {

	private static final Parallelism oneCpu = Parallelism.makeCpu(1);
	private static final Parallelism allCpus = Parallelism.makeCpu(Parallelism.getMaxNumCPUs());

	private static final double freeEnergyEpsilon = 1e-6;

	static {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();
	}

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

	/**
	 * A small peptide inhibitor affinity design
	 */
	public static MultiStateConfSpace affinity_6ov7_1mut6flex() {
		var confSpaces = new TestConfSpace.AffinityCompiled(
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.small.complex.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.small.design.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.small.target.ccsx"))
		);
		return new MultiStateConfSpace
			.Builder("complex", confSpaces.complex)
			.addMutableState("design", confSpaces.chainA)
			.addUnmutableState("target", confSpaces.chainB)
			.build();
	}

	private static void withPseudoCluster(int numMembers, Consumer<Cluster> block) {
		var exceptions = ClusterMember.launchPseudoCluster(numMembers, block);
		if (!exceptions.isEmpty()) {
			fail("Cluster threads encountered exceptions");
		}
	}

	private static Coffee makeCoffee(MultiStateConfSpace confSpace, Cluster cluster, Parallelism parallelism, long bytes) {

		// looser bounds makes the zmat calculation much faster, and the overall computation faster for small conf spaces
		var posInterDist = PosInterDist.DesmetEtAl1992;
		//var posInterDist = PosInterDist.TighterBounds;

		return new Coffee.Builder(confSpace)
			.setCluster(cluster)
			.setParallelism(parallelism)
			.setNodeDBMem(bytes)
			.configEachState(config -> {
				var stateConfSpace = (ConfSpace)config.state.confSpace;
				if (parallelism.numGpus > 0) {
					config.ecalc = new CudaConfEnergyCalculator(stateConfSpace, Structs.Precision.Float64, parallelism);
				} else {
					config.ecalc = new CPUConfEnergyCalculator(stateConfSpace);
				}
				config.posInterGen = new PosInterGen(posInterDist, null);
			})
			.build();
	}

	private static void bruteForceAll(Coffee coffee) {
		for (var seq : coffee.confSpace.seqSpace.getSequences()) {
			bruteForceFreeEnergies(coffee, seq);
		}
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
				var bigMath = new BigMath(coffee.mathContext).set(0);
				var progress = new Progress(rcs.getNumConformations().longValueExact());

				// batch confs together for speed
				var batchSize = stateInfo.config.ecalc.maxBatchSize();
				Consumer<List<int[]>> submitBatch = confBatch -> {
					tasks.submit(
						() -> stateInfo.zPaths(confBatch),
						zs -> {
							for (var z : zs) {
								bigMath.add(z);
							}
							progress.incrementProgress(zs.size());
						}
					);
				};

				// fill the batches and process them
				var confBatch = new ArrayList<int[]>(batchSize);
				for (var confList : allConfs) {
					confBatch.add(confList.stream()
						.mapToInt(i -> i)
						.toArray());
					if (confBatch.size() == batchSize) {
						submitBatch.accept(confBatch);
						confBatch = new ArrayList<>();
					}
				}
				if (confBatch.size() > 0) {
					submitBatch.accept(confBatch);
				}
				tasks.waitForFinish();

				// finally, calculate the free energy
				double g = gcalc.calc(bigMath.get());
				log("\tstate %s = %f", state.name, g);

				// HACKHACK: release the gpu streams if needed
				if (stateInfo.config.ecalc instanceof CudaConfEnergyCalculator) {
					((CudaConfEnergyCalculator)stateInfo.config.ecalc).recycleStreams();
				}
			}
		}
	}

	private void seqFreeEnergy(MultiStateConfSpace confSpace, Function<SeqSpace,Sequence> seqFunc, double[] freeEnergies, long bytes, double precision, int numMembers, Parallelism parallelism) {
		withPseudoCluster(numMembers, cluster -> {

			// get the sequence
			Coffee coffee = makeCoffee(confSpace, cluster, parallelism, bytes);
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

	// single-sequence test cases
	private void seqFreeEnergy_affinity_6ov7_1mut2flex_wt(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			new double[] { -1377.127950, -144.199934, -1187.667391 },
			bytes, precision, numMembers, parallelism
		);
	}
	private void seqFreeEnergy_affinity_6ov7_1mut2flex_ala(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ALA"),
			new double[] { -1363.561940, -132.431356, -1187.667391 },
			bytes, precision, numMembers, parallelism
		);
	}
	private void seqFreeEnergy_affinity_6ov7_1mut2flex_asn(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ASN"),
			new double[] { -1375.773406, -143.920583, -1187.667391 },
			bytes, precision, numMembers, parallelism
		);
	}
	private static void bruteForce_affinity_6ov7_1mut2flex(Parallelism parallelism) {
		bruteForceAll(makeCoffee(TestCoffee.affinity_6ov7_1mut2flex(), null, parallelism, 0));
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

	private void seqFreeEnergy_affinity_6ov7_1mut6flex_wt(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut6flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			new double[] { -1380.512795, -145.154178, -1190.413094 },
			bytes, precision, numMembers, parallelism
		);
	}
	private void seqFreeEnergy_affinity_6ov7_1mut6flex_ala(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut6flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ALA"),
			new double[] { -1366.686478, -133.531642, -1190.413094 },
			bytes, precision, numMembers, parallelism
		);
	}
	private void seqFreeEnergy_affinity_6ov7_1mut6flex_asn(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut6flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ASN"),
			new double[] { -1378.915230, -145.208459, -1190.413094 },
			bytes, precision, numMembers, parallelism
		);
	}
	private static void bruteForce_affinity_6ov7_1mut6flex(Parallelism parallelism) {
		bruteForceAll(makeCoffee(TestCoffee.affinity_6ov7_1mut6flex(), null, parallelism, 0));
		//sequence [6 GLN=gln]
		//	state complex = -1380.512795
		//	state design = -145.154178
		//	state target = -1190.413094
		//sequence [6 GLN=ALA]
		//	state complex = -1366.686478
		//	state design = -133.531642
		//	state target = -1190.413094
		//sequence [6 GLN=ASN]
		//	state complex = -1378.915230
		//	state design = -145.208459
		//	state target = -1190.413094
	}

	public static void main(String[] args) {
		var allGpus = Parallelism.make(3*4, 4, 1 /* ignored */);
		//bruteForce_affinity_6ov7_1mut2flex(allCpus);
		bruteForce_affinity_6ov7_1mut6flex(allGpus);
	}

	// TINY CONF SPACE

	// the basic test
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, 1, oneCpu);
	}

	// vary precision
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_001_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.01, 1, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_0001_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.001, 1, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_0_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.0001, 1, oneCpu);
	}

	// vary sequence
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_ala_01_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_ala(1024*1024, 0.1, 1, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_asn_01_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_asn(1024*1024, 0.1, 1, oneCpu);
	}

	// vary cluster members/parallelism
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_2x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, 2, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_2x2_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, 2, Parallelism.makeCpu(2));
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x2_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, 1, Parallelism.makeCpu(2));
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x4_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, 1, Parallelism.makeCpu(4));
	}

	// vary the memory
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x1_256k() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(256*1024, 0, 1, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x1_96k() { // min possible bytes, still enough space for perfect precision
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(96*1024, 0, 1, oneCpu);
	}


	// SMALL CONF SPACE

	// the basic test
	@Test public void seqFreeEnergy_affinity_6ov7_1mut6flex_wt_01_1x1_1m() {
		// TEMP
		//seqFreeEnergy_affinity_6ov7_1mut6flex_wt(1024*1024, 0.1, 1, oneCpu);
		seqFreeEnergy_affinity_6ov7_1mut6flex_wt(96*1024, 0.1, 1, oneCpu);
	}
}
