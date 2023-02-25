package edu.duke.cs.osprey.coffee;

import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.*;
import static edu.duke.cs.osprey.TestBase.*;
import static edu.duke.cs.osprey.tools.Log.log;
import static org.junit.jupiter.api.Assertions.fail;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.coffee.directors.*;
import edu.duke.cs.osprey.coffee.seqdb.SeqFreeEnergies;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.energy.compiled.*;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.*;
import edu.duke.cs.osprey.tools.MathTools.*;
import org.junit.jupiter.api.Test;

import java.math.MathContext;
import java.math.RoundingMode;
import java.time.Duration;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestCoffee {

	private static final Parallelism oneCpu = Parallelism.makeCpu(1);
	private static final Parallelism allCpus = Parallelism.makeCpu(Parallelism.getMaxNumCPUs());

	private static final double freeEnergyEpsilon = 1e-6;

	static {
		Cluster.fixHazelcastLogging();
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

	/**
	 * A medium peptide inhibitor affinity design
	 */
	public static MultiStateConfSpace affinity_6ov7_1mut11flex() {
		var confSpaces = new TestConfSpace.AffinityCompiled(
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.medium.complex.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.medium.design.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.medium.target.ccsx"))
		);
		return new MultiStateConfSpace
			.Builder("complex", confSpaces.complex)
			.addMutableState("design", confSpaces.chainA)
			.addUnmutableState("target", confSpaces.chainB)
			.build();
	}

	/**
	 * A small peptide inhibitor effinity design with ~50 sequences.
	 */
	public static MultiStateConfSpace affinity_6ov7_2mut4flex() {
		var confSpaces = new TestConfSpace.AffinityCompiled(
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.2m4f.complex.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.2m4f.design.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.2m4f.target.ccsx"))
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

	private static Coffee makeCoffee(MultiStateConfSpace confSpace, PosInterDist posInterDist, Double triplesThreshold, Cluster cluster, Parallelism parallelism, long bytes) {
		return new Coffee.Builder(confSpace)
			.setCluster(cluster)
			.setParallelism(parallelism)
			.setNodeDBMem(bytes)
			.setTripleCorrectionThreshold(triplesThreshold)
			.configEachState(config -> {
				config.posInterGen = new PosInterGen(posInterDist, null);
			})
			.build();
	}

	private static void bruteForceAll(Coffee coffee, BiConsumer<Coffee,Sequence> func) {
		for (var seq : coffee.confSpace.seqSpace.getSequences()) {
			func.accept(coffee, seq);
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

				// make the energy calculator
				try (var ecalc = ConfEnergyCalculator.makeBest(stateInfo.config.confSpace, coffee.parallelism)) {
					var batchSize = ecalc.maxBatchSize();

					// batch confs together for speed
					Consumer<List<ConfEnergyCalculator.MinimizationJob>> submitBatch = jobs -> {
						tasks.submit(
							() -> {
								ecalc.minimizeEnergies(jobs);
								return 42;
							},
							answer -> {
								for (var job : jobs) {
									bigMath.add(coffee.bcalc.calcPrecise(job.energy));
								}
								progress.incrementProgress(jobs.size());
							}
						);
					};

					// fill the batches and process them
					var jobs = new ArrayList<ConfEnergyCalculator.MinimizationJob>(batchSize);
					for (var confList : allConfs) {
						int[] conf = confList.stream()
							.mapToInt(i -> i)
							.toArray();
						var inters = stateInfo.config.posInterGen.all(stateInfo.config.confSpace, conf);
						jobs.add(new ConfEnergyCalculator.MinimizationJob(conf, inters));
						if (jobs.size() == batchSize) {
							submitBatch.accept(jobs);
							jobs = new ArrayList<>();
						}
					}
					if (jobs.size() > 0) {
						submitBatch.accept(jobs);
					}
					tasks.waitForFinish();
				}

				// finally, calculate the free energy
				double g = gcalc.calc(bigMath.get());
				log("\tstate %s = %f", state.name, g);
			}
		}
	}

	private static void gradientDescentFreeEnergies(Coffee coffee, Sequence seq) {

		log("sequence [%s]", seq);

		var gcalc = new FreeEnergyCalculator();

		try (var tasks = new ThreadPoolTaskExecutor()) {
			tasks.start(coffee.parallelism.numThreads);

			for (var state : coffee.confSpace.states) {
				var stateInfo = coffee.infos[state.index];

				try (var ctx = tasks.contextGroup()) {

					//var posInterDist = PosInterDist.DesmetEtAl1992;
					var posInterDist = PosInterDist.TighterBounds;

					// compute an emat
					var ecalc = new CPUConfEnergyCalculator(stateInfo.config.confSpace);
					var emat = new EmatCalculator.Builder(ecalc)
						.setPosInterDist(posInterDist)
						.setIncludeStaticStatic(true)
						.build()
						.calc();

					// compute the pfunc
					var confEcalc = new ConfEnergyCalculatorAdapter.Builder(ecalc, tasks)
						.setPosInterDist(posInterDist)
						.setIncludeStaticStatic(true)
						.build();
					var rcs = seq.makeRCs(stateInfo.config.confSpace);
					var pfunc = new GradientDescentPfunc(
						confEcalc,
						new ConfAStarTree.Builder(emat, rcs)
							.setTraditional()
							.build(),
						new ConfAStarTree.Builder(emat, rcs)
							.setTraditional()
							.build(),
						rcs.getNumConformations()
					);
					pfunc.init(0.01);
					pfunc.setInstanceId(0);
					pfunc.putTaskContexts(ctx);
					//pfunc.setReportProgress(true);
					pfunc.compute();

					// finally, calculate the free energy
					var values = pfunc.getValues();
					var g = gcalc.calc(new MathTools.BigDecimalBounds(
						values.calcLowerBound(),
						values.calcUpperBound()
					));
					log("\tstate %s = [%f,%f]", state.name, g.lower, g.upper);
				}
			}
		}
	}


	private void seqFreeEnergy(MultiStateConfSpace confSpace, Function<SeqSpace,Sequence> seqFunc, PosInterDist posInterDist, Double triplesThreshold, double[] freeEnergies, long bytes, double precision, int numMembers, Parallelism parallelism) {
		withPseudoCluster(numMembers, cluster -> {

			// get the sequence
			Coffee coffee = makeCoffee(confSpace, posInterDist, triplesThreshold, cluster, parallelism, bytes);
			Sequence seq = seqFunc.apply(coffee.confSpace.seqSpace);

			// run COFFEE
			var director = new SequenceDirector(coffee.confSpace, seq, precision, true);
			coffee.run(director);

			// check the free energies
			if (cluster.nodeId == 0) {
				for (var state : confSpace.states) {
					var g = director.getFreeEnergy(state);
					assertThat(state.name, g, isAbsoluteBound(freeEnergies[state.index], freeEnergyEpsilon));
					assertThat(state.name, g.size(), lessThanOrEqualTo(precision));
				}
			}
		});
	}

	private void seqFreeEnergy(MultiStateConfSpace confSpace, Function<SeqSpace,Sequence> seqFunc, PosInterDist posInterDist, Double triplesThreshold, DoubleBounds[] freeEnergies, long bytes, double precision, int numMembers, Parallelism parallelism) {
		withPseudoCluster(numMembers, cluster -> {

			// get the sequence
			Coffee coffee = makeCoffee(confSpace, posInterDist, triplesThreshold, cluster, parallelism, bytes);
			Sequence seq = seqFunc.apply(coffee.confSpace.seqSpace);

			// run COFFEE
			var director = new SequenceDirector(coffee.confSpace, seq, precision, true);
			coffee.run(director);

			// check the free energies
			if (cluster.nodeId == 0) {
				for (var state : confSpace.states) {
					var g = director.getFreeEnergy(state);
					assertThat(state.name, g, intersectsAbsolutely(freeEnergies[state.index], freeEnergyEpsilon));
					assertThat(state.name, g.size(), lessThanOrEqualTo(precision));
				}
			}
		});
	}

	// single-sequence test cases
	private void seqFreeEnergy_affinity_6ov7_1mut2flex_wt(long bytes, double precision, boolean triples, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			PosInterDist.DesmetEtAl1992,
			triples ? 10.0 : null,
			new double[] { -1377.127950, -144.199934, -1187.667391 },
			bytes, precision, numMembers, parallelism
		);
	}
	private void seqFreeEnergy_affinity_6ov7_1mut2flex_ala(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ALA"),
			PosInterDist.DesmetEtAl1992,
			null,
			new double[] { -1363.561940, -132.431356, -1187.667391 },
			bytes, precision, numMembers, parallelism
		);
	}
	private void seqFreeEnergy_affinity_6ov7_1mut2flex_asn(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut2flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ASN"),
			PosInterDist.DesmetEtAl1992,
			null,
			new double[] { -1375.773406, -143.920583, -1187.667391 },
			bytes, precision, numMembers, parallelism
		);
	}
	private static void bruteForce_affinity_6ov7_1mut2flex(Parallelism parallelism) {
		bruteForceAll(
			makeCoffee(TestCoffee.affinity_6ov7_1mut2flex(), PosInterDist.DesmetEtAl1992, null, null, parallelism, 0),
			TestCoffee::bruteForceFreeEnergies
		);
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

	private void seqFreeEnergy_affinity_6ov7_1mut6flex_wt(long bytes, double precision, boolean triples, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut6flex(),
			seqSpace -> seqSpace.makeWildTypeSequence(),
			PosInterDist.TighterBounds, // tighter bounds here cuts the runtime in half!
			triples ? 1.0 : null,
			new DoubleBounds[] {
				new DoubleBounds(-1380.529262, -1380.523417),
				new DoubleBounds(-145.154147, -145.153834),
				new DoubleBounds(-1190.451347, -1190.446539)
			},
			bytes, precision, numMembers, parallelism
		);
	}
	private void seqFreeEnergy_affinity_6ov7_1mut6flex_ala(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut6flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ALA"),
			PosInterDist.TighterBounds,
			null,
			new DoubleBounds[] {
				new DoubleBounds(-1366.733962, -1366.728450),
				new DoubleBounds(-133.531642, -133.531642),
				new DoubleBounds(-1190.451347, -1190.446539)
			},
			bytes, precision, numMembers, parallelism
		);
	}
	private void seqFreeEnergy_affinity_6ov7_1mut6flex_asn(long bytes, double precision, int numMembers, Parallelism parallelism) {
		seqFreeEnergy(
			TestCoffee.affinity_6ov7_1mut6flex(),
			seqSpace -> seqSpace.makeWildTypeSequence().set("6 GLN", "ASN"),
			PosInterDist.TighterBounds,
			null,
			new DoubleBounds[] {
				new DoubleBounds(-1378.941314, -1378.935465),
				new DoubleBounds(-145.208390, -145.208390),
				new DoubleBounds(-1190.451347, -1190.446539)
			},
			bytes, precision, numMembers, parallelism
		);
	}
	private static void bruteForce_affinity_6ov7_1mut6flex(Parallelism parallelism) {
		bruteForceAll(
			makeCoffee(TestCoffee.affinity_6ov7_1mut6flex(), PosInterDist.TighterBounds, null, null, parallelism, 0),
			TestCoffee::gradientDescentFreeEnergies
		);
		//sequence [6 GLN=gln]
		//	state complex = [-1380.529262,-1380.523417]
		//	state design = [-145.154147,-145.153834]
		//	state target = [-1190.451347,-1190.446539]
		//sequence [6 GLN=ALA]
		//	state complex = [-1366.733962,-1366.728450]
		//	state design = [-133.531642,-133.531642]
		//	state target = [-1190.451347,-1190.446539]

		//sequence [6 GLN=ASN]
		//	state complex = [-1378.941314,-1378.935465]
		//	state design = [-145.208390,-145.208390]
		//	state target = [-1190.451347,-1190.446539]
	}

	public static void main(String[] args) {
		bruteForce_affinity_6ov7_1mut2flex(allCpus);
		bruteForce_affinity_6ov7_1mut6flex(allCpus);
		design_affinity_6ov7_2mut4flex_precise();
	}

	// TINY CONF SPACE

	// the basic test
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, false, 1, oneCpu);
	}

	// vary precision
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_001_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.01, false, 1, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_0001_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.001, false, 1, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_0_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.0001, false, 1, oneCpu);
	}

	// vary sequence
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_ala_01_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_ala(1024*1024, 0.1, 1, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_asn_01_1x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_asn(1024*1024, 0.1, 1, oneCpu);
	}

	// vary cluster members/parallelism
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x2_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, false, 1, Parallelism.makeCpu(2));
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x4_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, false, 1, Parallelism.makeCpu(4));
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_2x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, false, 2, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_2x2_1m() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, false, 2, Parallelism.makeCpu(2));
	}

	// vary the memory
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_0_1x1_256k() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(256*1024, 0, false, 1, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_0_1x1_96k() { // min possible bytes, still enough space for perfect precision
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(96*1024, 0, false, 1, oneCpu);
	}

	// add triples
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x1_1m_triples() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, true, 1, oneCpu);
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_1x4_1m_triples() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, true, 1, Parallelism.makeCpu(4));
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut2flex_wt_01_2x1_1m_triples() {
		seqFreeEnergy_affinity_6ov7_1mut2flex_wt(1024*1024, 0.1, true, 2, oneCpu);
	}


	// SMALL CONF SPACE

	// the basic test
	@Test public void seqFreeEnergy_affinity_6ov7_1mut6flex_wt_01_1x4_1m() {
		seqFreeEnergy_affinity_6ov7_1mut6flex_wt(1024*1024, 1.0, false, 1, allCpus);
	}

	// vary the cluster members/parallelism
	@Test public void seqFreeEnergy_affinity_6ov7_1mut6flex_wt_01_2x2_1m() {
		seqFreeEnergy_affinity_6ov7_1mut6flex_wt(1024*1024, 1.0, false, 2, Parallelism.makeCpu(2));
	}
	@Test public void seqFreeEnergy_affinity_6ov7_1mut6flex_wt_01_4x1_1m() {
		seqFreeEnergy_affinity_6ov7_1mut6flex_wt(1024*1024, 1.0, false, 4, oneCpu);
	}

	// add triples
	//@Test // this test takes 12 minutes on my laptop... waay too long to leave on by default
	public void seqFreeEnergy_affinity_6ov7_1mut6flex_wt_01_1x4_1m_triples() {
		seqFreeEnergy_affinity_6ov7_1mut6flex_wt(1024*1024, 1.0, true, 1, allCpus);
	}
	
	private void assertFreeEnergies_2m4f_complex(Collection<SeqFreeEnergies> obsgs, int complexi) {
		
		var expgs = new HashMap<String,Double>();

		// expected values computed by running COFFEE on each sequence separately to gWidthMax = 1e-3
		// and with a SeqDB precision of 4096, to accurately compute the high-energy sequences
		expgs.put("8 THR=thr 10 VAL=val", -1378.275);

		expgs.put("8 THR=ASN 10 VAL=val", -1379.407);
		expgs.put("8 THR=ASP 10 VAL=val", -1381.103);
		expgs.put("8 THR=GLN 10 VAL=val", -1383.748);
		expgs.put("8 THR=GLU 10 VAL=val", -1381.936);
		expgs.put("8 THR=GLY 10 VAL=val", -1369.435);
		expgs.put("8 THR=SER 10 VAL=val", -1376.251);

		expgs.put("8 THR=thr 10 VAL=ALA", -1371.377);
		expgs.put("8 THR=thr 10 VAL=ILE", -1370.219);
		expgs.put("8 THR=thr 10 VAL=LEU", -1348.772);
		expgs.put("8 THR=thr 10 VAL=PHE", -702.158);
		expgs.put("8 THR=thr 10 VAL=TRP", 2132.812);
		expgs.put("8 THR=thr 10 VAL=TYR", 1513.490);

		expgs.put("8 THR=ASN 10 VAL=ALA", -1372.557);
		expgs.put("8 THR=ASN 10 VAL=ILE", -1371.100);
		expgs.put("8 THR=ASN 10 VAL=LEU", -1353.375);
		expgs.put("8 THR=ASN 10 VAL=PHE", -706.786);
		expgs.put("8 THR=ASN 10 VAL=TRP", 2129.810);
		expgs.put("8 THR=ASN 10 VAL=TYR", 1504.753);

		expgs.put("8 THR=ASP 10 VAL=ALA", -1374.422);
		expgs.put("8 THR=ASP 10 VAL=ILE", -1372.843);
		expgs.put("8 THR=ASP 10 VAL=LEU", -1353.972);
		expgs.put("8 THR=ASP 10 VAL=PHE", -706.888);
		expgs.put("8 THR=ASP 10 VAL=TRP", 2125.895);
		expgs.put("8 THR=ASP 10 VAL=TYR", 1504.968);

		expgs.put("8 THR=GLN 10 VAL=ALA", -1376.819);
		expgs.put("8 THR=GLN 10 VAL=ILE", -1375.701);
		expgs.put("8 THR=GLN 10 VAL=LEU", -1356.205);
		expgs.put("8 THR=GLN 10 VAL=PHE", -709.313);
		expgs.put("8 THR=GLN 10 VAL=TRP", 2124.343);
		expgs.put("8 THR=GLN 10 VAL=TYR", 1502.252);

		expgs.put("8 THR=GLU 10 VAL=ALA", -1375.200);
		expgs.put("8 THR=GLU 10 VAL=ILE", -1373.810);
		expgs.put("8 THR=GLU 10 VAL=LEU", -1355.209);
		expgs.put("8 THR=GLU 10 VAL=PHE", -708.984);
		expgs.put("8 THR=GLU 10 VAL=TRP", 2124.959);
		expgs.put("8 THR=GLU 10 VAL=TYR", 1502.685);

		expgs.put("8 THR=GLY 10 VAL=ALA", -1362.851);
		expgs.put("8 THR=GLY 10 VAL=ILE", -1361.027);
		expgs.put("8 THR=GLY 10 VAL=LEU", -1342.512);
		expgs.put("8 THR=GLY 10 VAL=PHE", -699.969);
		expgs.put("8 THR=GLY 10 VAL=TRP", 2137.699);
		expgs.put("8 THR=GLY 10 VAL=TYR", 1510.260);

		expgs.put("8 THR=SER 10 VAL=ALA", -1369.563);
		expgs.put("8 THR=SER 10 VAL=ILE", -1368.029);
		expgs.put("8 THR=SER 10 VAL=LEU", -1349.555);
		expgs.put("8 THR=SER 10 VAL=PHE", -704.956);
		expgs.put("8 THR=SER 10 VAL=TRP", 2130.538);
		expgs.put("8 THR=SER 10 VAL=TYR", 1506.171);

		final double epsilon = 1e-3;

		for (var seqg : obsgs) {
			double expg = expgs.get(seqg.seq.toString());
			assertThat("unexpected sequence: " + seqg.seq, expg, is(not(nullValue())));
			assertThat("wrong free energy for: " + seqg.seq, seqg.freeEnergies[complexi], isAbsoluteBound(expg, epsilon));
		}
	}

	private void assertFreeEnergies_2m4f_design(Collection<SeqFreeEnergies> obsgs, int designi) {

		var expgs = new HashMap<String,Double>();

		expgs.put("8 THR=thr 10 VAL=val", -143.431);

		expgs.put("8 THR=ASN 10 VAL=val", -148.725);
		expgs.put("8 THR=ASP 10 VAL=val", -147.593);
		expgs.put("8 THR=GLN 10 VAL=val", -149.007);
		expgs.put("8 THR=GLU 10 VAL=val", -147.765);
		expgs.put("8 THR=GLY 10 VAL=val", -137.280);
		expgs.put("8 THR=SER 10 VAL=val", -142.320);

		expgs.put("8 THR=thr 10 VAL=ALA", -140.970);
		expgs.put("8 THR=thr 10 VAL=ILE", -141.569);
		expgs.put("8 THR=thr 10 VAL=LEU", -143.159);
		expgs.put("8 THR=thr 10 VAL=PHE", -142.346);
		expgs.put("8 THR=thr 10 VAL=TRP", -147.089);
		expgs.put("8 THR=thr 10 VAL=TYR", -146.971);

		expgs.put("8 THR=ASN 10 VAL=ALA", -146.348);
		expgs.put("8 THR=ASN 10 VAL=ILE", -146.877);
		expgs.put("8 THR=ASN 10 VAL=LEU", -148.286);
		expgs.put("8 THR=ASN 10 VAL=PHE", -148.076);
		expgs.put("8 THR=ASN 10 VAL=TRP", -152.489);
		expgs.put("8 THR=ASN 10 VAL=TYR", -152.787);

		expgs.put("8 THR=ASP 10 VAL=ALA", -145.224);
		expgs.put("8 THR=ASP 10 VAL=ILE", -145.740);
		expgs.put("8 THR=ASP 10 VAL=LEU", -147.204);
		expgs.put("8 THR=ASP 10 VAL=PHE", -146.404);
		expgs.put("8 THR=ASP 10 VAL=TRP", -151.039);
		expgs.put("8 THR=ASP 10 VAL=TYR", -151.053);

		expgs.put("8 THR=GLN 10 VAL=ALA", -146.489);
		expgs.put("8 THR=GLN 10 VAL=ILE", -147.141);
		expgs.put("8 THR=GLN 10 VAL=LEU", -148.775);
		expgs.put("8 THR=GLN 10 VAL=PHE", -148.251);
		expgs.put("8 THR=GLN 10 VAL=TRP", -152.825);
		expgs.put("8 THR=GLN 10 VAL=TYR", -152.965);

		expgs.put("8 THR=GLU 10 VAL=ALA", -145.423);
		expgs.put("8 THR=GLU 10 VAL=ILE", -145.906);
		expgs.put("8 THR=GLU 10 VAL=LEU", -147.484);
		expgs.put("8 THR=GLU 10 VAL=PHE", -147.278);
		expgs.put("8 THR=GLU 10 VAL=TRP", -152.043);
		expgs.put("8 THR=GLU 10 VAL=TYR", -152.471);

		expgs.put("8 THR=GLY 10 VAL=ALA", -135.107);
		expgs.put("8 THR=GLY 10 VAL=ILE", -135.442);
		expgs.put("8 THR=GLY 10 VAL=LEU", -136.811);
		expgs.put("8 THR=GLY 10 VAL=PHE", -136.723);
		expgs.put("8 THR=GLY 10 VAL=TRP", -140.940);
		expgs.put("8 THR=GLY 10 VAL=TYR", -141.269);

		expgs.put("8 THR=SER 10 VAL=ALA", -140.033);
		expgs.put("8 THR=SER 10 VAL=ILE", -140.475);
		expgs.put("8 THR=SER 10 VAL=LEU", -141.909);
		expgs.put("8 THR=SER 10 VAL=PHE", -141.767);
		expgs.put("8 THR=SER 10 VAL=TRP", -146.227);
		expgs.put("8 THR=SER 10 VAL=TYR", -146.424);

		final double epsilon = 1e-3;

		for (var seqg : obsgs) {
			double expg = expgs.get(seqg.seq.toString());
			assertThat("unexpected sequence: " + seqg.seq, expg, is(not(nullValue())));
			assertThat("wrong free energy for: " + seqg.seq, seqg.freeEnergies[designi], isAbsoluteBound(expg, epsilon));
		}
	}

	private static void design_affinity_6ov7_2mut4flex_precise() {

		var confSpace = TestCoffee.affinity_6ov7_2mut4flex();

		var coffee = new Coffee.Builder(confSpace)
			.setParallelism(allCpus)
			.setSeqDBMathContext(new MathContext(4096, RoundingMode.HALF_UP))
			.configEachState(config -> {
				config.posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
			})
			.build();

		var seqs = Arrays.asList(
			confSpace.seqSpace.makeSequence("thr", "val"),
			confSpace.seqSpace.makeSequence("ASN", "val"),
			confSpace.seqSpace.makeSequence("ASP", "val"),
			confSpace.seqSpace.makeSequence("GLN", "val"),
			confSpace.seqSpace.makeSequence("GLU", "val"),
			confSpace.seqSpace.makeSequence("GLY", "val"),
			confSpace.seqSpace.makeSequence("SER", "val"),
			confSpace.seqSpace.makeSequence("thr", "ALA"),
			confSpace.seqSpace.makeSequence("thr", "ILE"),
			confSpace.seqSpace.makeSequence("thr", "LEU"),
			confSpace.seqSpace.makeSequence("thr", "PHE"),
			confSpace.seqSpace.makeSequence("thr", "TRP"),
			confSpace.seqSpace.makeSequence("thr", "TYR"),
			confSpace.seqSpace.makeSequence("ASN", "ALA"),
			confSpace.seqSpace.makeSequence("ASN", "ILE"),
			confSpace.seqSpace.makeSequence("ASN", "LEU"),
			confSpace.seqSpace.makeSequence("ASN", "PHE"),
			confSpace.seqSpace.makeSequence("ASN", "TRP"),
			confSpace.seqSpace.makeSequence("ASN", "TYR"),
			confSpace.seqSpace.makeSequence("ASP", "ALA"),
			confSpace.seqSpace.makeSequence("ASP", "ILE"),
			confSpace.seqSpace.makeSequence("ASP", "LEU"),
			confSpace.seqSpace.makeSequence("ASP", "PHE"),
			confSpace.seqSpace.makeSequence("ASP", "TRP"),
			confSpace.seqSpace.makeSequence("ASP", "TYR"),
			confSpace.seqSpace.makeSequence("GLN", "ALA"),
			confSpace.seqSpace.makeSequence("GLN", "ILE"),
			confSpace.seqSpace.makeSequence("GLN", "LEU"),
			confSpace.seqSpace.makeSequence("GLN", "PHE"),
			confSpace.seqSpace.makeSequence("GLN", "TRP"),
			confSpace.seqSpace.makeSequence("GLN", "TYR"),
			confSpace.seqSpace.makeSequence("GLU", "ALA"),
			confSpace.seqSpace.makeSequence("GLU", "ILE"),
			confSpace.seqSpace.makeSequence("GLU", "LEU"),
			confSpace.seqSpace.makeSequence("GLU", "PHE"),
			confSpace.seqSpace.makeSequence("GLU", "TRP"),
			confSpace.seqSpace.makeSequence("GLU", "TYR"),
			confSpace.seqSpace.makeSequence("GLY", "ALA"),
			confSpace.seqSpace.makeSequence("GLY", "ILE"),
			confSpace.seqSpace.makeSequence("GLY", "LEU"),
			confSpace.seqSpace.makeSequence("GLY", "PHE"),
			confSpace.seqSpace.makeSequence("GLY", "TRP"),
			confSpace.seqSpace.makeSequence("GLY", "TYR"),
			confSpace.seqSpace.makeSequence("SER", "ALA"),
			confSpace.seqSpace.makeSequence("SER", "ILE"),
			confSpace.seqSpace.makeSequence("SER", "LEU"),
			confSpace.seqSpace.makeSequence("SER", "PHE"),
			confSpace.seqSpace.makeSequence("SER", "TRP"),
			confSpace.seqSpace.makeSequence("SER", "TYR")
		);

		for (var state : confSpace.states) {
			log("state: %s", state.name);

			var stateSeqs = Collections.singletonList((Sequence)null);
			if (state.isSequenced) {
				stateSeqs = seqs;
			}

			for (var seq : stateSeqs) {

				var director = new PfuncDirector.Builder(confSpace, state)
					.setSequence(seq)
					.setTiming(Timing.Precise)
					.setGWidthMax(1e-3)
					.build();
				coffee.run(director);

				var g = director.getFreeEnergy();
				assertThat(g.size(), lessThanOrEqualTo(1e-3));

				log("\texpgs.put(\"%s\", %.3f);", seq, g.lower);
			}
		}
	}

	private void design_affinity_6ov7_2mut4flex(int numMembers, int numThreads) {
		withPseudoCluster(numMembers, cluster -> {
			var confSpace = TestCoffee.affinity_6ov7_2mut4flex();

			Coffee coffee = new Coffee.Builder(confSpace)
				.setCluster(cluster)
				.setParallelism(Parallelism.makeCpu(numThreads))
				.configEachState(config -> {
					config.posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
				})
				.build();

			var director = new AffinityDirector.Builder(confSpace, "complex", "design", "target")
				.setK(28) // just the sequences that are quick to get... some of the other sequences have high energies and VERY loose bounds
				.setMaxSimultaneousMutations(null)
				.setTiming(Timing.Precise)
				.build();
			coffee.run(director);

			if (cluster.nodeId == 0) {
				assertThat(director.bestSeqs.size(), is(28));
				assertFreeEnergies_2m4f_complex(director.bestSeqs, confSpace.getState("complex").index);
				assertFreeEnergies_2m4f_design(director.bestSeqs, confSpace.getState("design").index);
				assertThat(director.targetFreeEnergy, isAbsoluteBound(-1187.609, 1e-3));
			}
		});
	}

	@Test public void design_affinity_6ov7_2mut4flex_1x4() { design_affinity_6ov7_2mut4flex(1, 4); }
	@Test public void design_affinity_6ov7_2mut4flex_2x2() { design_affinity_6ov7_2mut4flex(2, 2); }
	@Test public void design_affinity_6ov7_2mut4flex_4x1() { design_affinity_6ov7_2mut4flex(4, 1); }

	@Test
	public void design_affinity_6ov7_2mut4flex_1x4_1mut() {
		var confSpace = TestCoffee.affinity_6ov7_2mut4flex();

		Coffee coffee = new Coffee.Builder(confSpace)
			.setParallelism(Parallelism.makeCpu(4))
			.configEachState(config -> {
				config.posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
			})
			.build();

		var director = new AffinityDirector.Builder(confSpace, "complex", "design", "target")
			.setK(11) // the other 2 single-mutation sequences have high energies and take forever to find
			.setMaxSimultaneousMutations(1)
			.setTiming(Timing.Precise)
			.build();
		coffee.run(director);

		assertThat(director.bestSeqs.size(), is(11));
		for (var seqg : director.bestSeqs) {
			assertThat(seqg.seq.toString(), seqg.seq.countMutations(), lessThanOrEqualTo(1));
		}

		assertFreeEnergies_2m4f_complex(director.bestSeqs, confSpace.getState("complex").index);
		assertFreeEnergies_2m4f_design(director.bestSeqs, confSpace.getState("design").index);
		assertThat(director.targetFreeEnergy, isAbsoluteBound(-1187.609, 1e-3));
	}


	@Test
	public void design_kstar_6ov7_2mut4flex_1x4_nothreshold() {
		var confSpace = TestCoffee.affinity_6ov7_2mut4flex();

		Coffee coffee = new Coffee.Builder(confSpace)
			.setParallelism(Parallelism.makeCpu(4))
			.configEachState(config -> {
				config.posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
			})
			.build();

		var director = new KStarDirector.Builder(confSpace, "complex", "design", "target")
			.setMaxSimultaneousMutations(null)
			.setStabilityThreshold(null)
			.setTiming(Timing.Precise)
			.setReportStateProgress(false)
			.build();

		// use just a few sequences for this test
		director.sequences.add(confSpace.seqSpace.makeSequence("GLN", "val"));
		director.sequences.add(confSpace.seqSpace.makeSequence("ASN", "TYR"));
		director.sequences.add(confSpace.seqSpace.makeSequence("GLU", "ILE"));

		coffee.run(director);

		assertThat(director.wildTypeG, is(nullValue()));
		assertThat(director.sequenceGs.size(), is(3));
		assertFreeEnergies_2m4f_complex(director.sequenceGs, confSpace.getState("complex").index);
		assertFreeEnergies_2m4f_design(director.sequenceGs, confSpace.getState("design").index);
		assertThat(director.targetG, isAbsoluteBound(-1187.609, 1e-3));
	}

	@Test
	public void design_kstar_6ov7_2mut4flex_1x4_threshold5() {
		var confSpace = TestCoffee.affinity_6ov7_2mut4flex();
		var complexi = confSpace.getState("complex").index;
		var designi = confSpace.getState("design").index;

		Coffee coffee = new Coffee.Builder(confSpace)
			.setParallelism(Parallelism.makeCpu(4))
			.configEachState(config -> {
				config.posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
			})
			.build();

		var director = new KStarDirector.Builder(confSpace, "complex", "design", "target")
			.setMaxSimultaneousMutations(null)
			.setStabilityThreshold(5.0)
			.setTiming(Timing.Precise)
			.setReportStateProgress(false)
			.build();

		// use just a few sequences for this test
		director.sequences.add(confSpace.seqSpace.makeSequence("GLN", "val"));
		director.sequences.add(confSpace.seqSpace.makeSequence("GLY", "ALA"));
		director.sequences.add(confSpace.seqSpace.makeSequence("GLU", "ILE"));

		coffee.run(director);

		assertThat(director.sequenceGs.size(), is(3));
		assertFreeEnergies_2m4f_complex(List.of(director.wildTypeG), complexi);
		assertFreeEnergies_2m4f_design(List.of(director.wildTypeG), designi);
		assertFreeEnergies_2m4f_design(director.sequenceGs, designi);
		assertThat(director.targetG, isAbsoluteBound(-1187.609, 1e-3));

		// GLY, ALA should fail the stability filter
		assertThat(director.sequenceGs.get(1).freeEnergies[complexi], is(nullValue()));

		// but not the others
		assertFreeEnergies_2m4f_complex(List.of(director.sequenceGs.get(0), director.sequenceGs.get(2)), complexi);
	}

	/**
	 * not a regular test, this design is waaaaay too big to test all the time
	 * but it's useful for doing in-depth debugging when you need bigger designs
	 */
	//@Test
	public void pfuncs_6ov7_2mut4flex_1x4() {

		var stateConfSpace = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.2m4f.complex.ccsx"));
		var confSpace = new MultiStateConfSpace.Builder("state", stateConfSpace)
			.build();
		var state = confSpace.states.get(0);

		Coffee coffee = new Coffee.Builder(confSpace)
			.setParallelism(Parallelism.makeCpu(4))
			.configEachState(config -> {
				config.posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
				//config.posInterGen = new PosInterGen(PosInterDist.TighterBounds, null);
			})
			.setNodeStatsReportingInterval(Duration.ofMillis(1_000))
			.build();

		var director = new PfuncsDirector.Builder(confSpace, state)
			//.setTiming(Timing.Precise)
			.setTiming(Timing.Efficient)
			.setReportProgress(true)
			//.addSequences(confSpace.seqSpace.getSequences(1))
			.addSequence(state.confSpace.seqSpace().makeSequence("THR", "TRP"))
			.build();

		coffee.run(director);
	}
}
