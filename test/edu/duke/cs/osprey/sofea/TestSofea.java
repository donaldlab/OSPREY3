package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.TestBase.isAbsoluteBound;
import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase.TempFile;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.lute.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.BigIntegerBounds;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.Streams;
import org.junit.Test;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.function.Consumer;


public class TestSofeaLute {

	// TODO: test parallel settings

	private static final Parallelism fullCPUParallelism = Parallelism.makeCpu(Parallelism.getMaxNumCPUs());
	private static final File tmpdir = new File(System.getProperty("java.io.tmpdir"), "testSofea");
	private static final BoltzmannCalculator bcalc = new BoltzmannCalculator(new MathContext(32, RoundingMode.HALF_UP));
	private static final double epsilonG = 1e-3;

	static {
		tmpdir.mkdirs();
	}


	@Test
	public void test_Binding1CC8Flex4_Standard_LeafCounts() {
		Design design = Designs.Binding1CC8Flex4_Standard.get();
		assertLeafCounts(design);
	}
	@Test
	public void test_Binding1CC8Flex4_Standard_ZValueBounds() {
		Design design = Designs.Binding1CC8Flex4_Standard.get();
		assertZValueBounds(design);
	}
	@Test
	public void test_Binding1CC8Flex4_Standard_ZSumBounds() {
		Design design = Designs.Binding1CC8Flex4_Standard.get();
		assertZSumBounds(design);
	}
	@Test
	public void test_Binding1CC8Flex4_Standard_CalcG() {
		Design design = Designs.Binding1CC8Flex4_Standard.get();
		assertGStates(design, Collections.emptyList(), -17.356, -79.225, -49.195);
	}


	@Test
	public void test_Stability1CC8Mut2Flex2_Standard_LeafCounts() {
		Design design = Designs.Stability1CC8Mut2Flex2_Standard.get();
		assertLeafCounts(design);
	}
	@Test
	public void test_Stability1CC8Mut2Flex2_Standard_CalcG() {
		Design design = Designs.Stability1CC8Mut2Flex2_Standard.get();
		assertGStates(design, Arrays.asList("SER", "GLY"), -17.356);
		assertGStates(design, Arrays.asList("VAL", "GLY"), -8.825);
		assertGStates(design, Arrays.asList("SER", "VAL"), -12.050);
		assertGStates(design, Arrays.asList("VAL", "VAL"), -3.617);
	}
	private static void assertResults_Stability1CC8Mut2Flex2_Standard(Results results) {
		results.assertGSequenced(Arrays.asList("SER", "GLY"), -17.356);
		results.assertGSequenced(Arrays.asList("VAL", "GLY"), -8.825);
		results.assertGSequenced(Arrays.asList("SER", "VAL"), -12.050);
		results.assertGSequenced(Arrays.asList("VAL", "VAL"), -3.617);
	}
	@Test
	public void test_Stability1CC8Mut2Flex2_Standard_SingleSweep() {
		sweepUntilExhaustion(
			Designs.Stability1CC8Mut2Flex2_Standard.get(),
			Double.POSITIVE_INFINITY,
			1024*1024,
			TestSofeaLute::assertResults_Stability1CC8Mut2Flex2_Standard
		);
	}
	@Test
	public void test_Stability1CC8Mut2Flex2_Standard_MultiSweepHiMem() {
		sweepUntilExhaustion(
			Designs.Stability1CC8Mut2Flex2_Standard.get(),
			50.0,
			1024*1024,
			TestSofeaLute::assertResults_Stability1CC8Mut2Flex2_Standard
		);
	}
	@Test
	public void test_Stability1CC8Mut2Flex2_Standard_MultiSweepLoMem() {
		sweepUntilExhaustion(
			Designs.Stability1CC8Mut2Flex2_Standard.get(),
			3.0,
			200,
			TestSofeaLute::assertResults_Stability1CC8Mut2Flex2_Standard
		);
	}


	@Test
	public void test_Binding1CC8Mut2Flex2_Standard_LeafCounts() {
		Design design = Designs.Binding1CC8Mut2Flex2_Standard.get();
		assertLeafCounts(design);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_Standard_ZValueBounds() {
		Design design = Designs.Binding1CC8Mut2Flex2_Standard.get();
		assertZValueBounds(design);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_Standard_ZSumBounds() {
		Design design = Designs.Binding1CC8Mut2Flex2_Standard.get();
		assertZSumBounds(design);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_Standard_CalcG() {
		Design design = Designs.Binding1CC8Mut2Flex2_Standard.get();
		assertGStates(design, Arrays.asList("SER", "GLY"), -17.356, -79.225, -49.195);
		assertGStates(design, Arrays.asList("VAL", "GLY"), -8.825, Double.POSITIVE_INFINITY, -49.195);
		assertGStates(design, Arrays.asList("SER", "VAL"), -12.050, Double.POSITIVE_INFINITY, -49.195);
		assertGStates(design, Arrays.asList("VAL", "VAL"), -3.617, Double.POSITIVE_INFINITY, -49.195);
	}
	private static void assertResults_Binding1CC8Mut2Flex2_Standard(Results results) {
		results.assertGUnsequenced(-49.195);
		results.assertGSequenced(Arrays.asList("SER", "GLY"), -17.356, -79.225);
		results.assertGSequenced(Arrays.asList("VAL", "GLY"), -8.825, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("SER", "VAL"), -12.050, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("VAL", "VAL"), -3.617, Double.POSITIVE_INFINITY);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_Standard_SingleSweep() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut2Flex2_Standard.get(),
			Double.POSITIVE_INFINITY,
			1024*1024,
			TestSofeaLute::assertResults_Binding1CC8Mut2Flex2_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_Standard_MultiSweepHiMem() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut2Flex2_Standard.get(),
			50.0,
			1024*1024,
			TestSofeaLute::assertResults_Binding1CC8Mut2Flex2_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_Standard_MultiSweepLoMem() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut2Flex2_Standard.get(),
			3.0,
			4*1024,
			TestSofeaLute::assertResults_Binding1CC8Mut2Flex2_Standard
		);
	}


	// too big for brute-force tests
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_CalcG() {
		Design design = Designs.Binding1CC8Mut3Flex4_Standard.get();
		assertGStates(design, Arrays.asList("SER", "GLY", "LYS"), -32.567, -148.509, -94.449);
		assertGStates(design, Arrays.asList("VAL", "GLY", "LYS"), -24.078, Double.POSITIVE_INFINITY, -94.449);
		assertGStates(design, Arrays.asList("SER", "VAL", "LYS"), -27.293, Double.POSITIVE_INFINITY, -94.449);
		assertGStates(design, Arrays.asList("VAL", "VAL", "LYS"), -18.904, Double.POSITIVE_INFINITY, -94.449);
		assertGStates(design, Arrays.asList("SER", "GLY", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, -94.449);
		assertGStates(design, Arrays.asList("VAL", "GLY", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, -94.449);
		assertGStates(design, Arrays.asList("SER", "VAL", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, -94.449);
		assertGStates(design, Arrays.asList("VAL", "VAL", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, -94.449);
	}
	private static void assertResults_Binding1CC8Mut3Flex4_Standard(Results results) {
		results.assertGUnsequenced(-94.449);
		results.assertGSequenced(Arrays.asList("SER", "GLY", "LYS"), -32.567, -148.509);
		results.assertGSequenced(Arrays.asList("VAL", "GLY", "LYS"), -24.078, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("SER", "VAL", "LYS"), -27.293, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("VAL", "VAL", "LYS"), -18.904, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("SER", "GLY", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("VAL", "GLY", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("SER", "VAL", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("VAL", "VAL", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
	}
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_SingleSweep() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut3Flex4_Standard.get(),
			Double.POSITIVE_INFINITY,
			1024*1024,
			TestSofeaLute::assertResults_Binding1CC8Mut3Flex4_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_SingleSweep_2Threads() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut3Flex4_Standard.get(),
			Double.POSITIVE_INFINITY,
			1024*1024,
			2,
			TestSofeaLute::assertResults_Binding1CC8Mut3Flex4_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_SingleSweep_4Threads() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut3Flex4_Standard.get(),
			Double.POSITIVE_INFINITY,
			1024*1024,
			4,
			TestSofeaLute::assertResults_Binding1CC8Mut3Flex4_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_MultiSweepHiMem() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut3Flex4_Standard.get(),
			50.0,
			1024*1024,
			TestSofeaLute::assertResults_Binding1CC8Mut3Flex4_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_MultiSweepHiMem_2Threads() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut3Flex4_Standard.get(),
			50.0,
			1024*1024,
			2,
			TestSofeaLute::assertResults_Binding1CC8Mut3Flex4_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_MultiSweepHiMem_4Threads() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut3Flex4_Standard.get(),
			50.0,
			1024*1024,
			4,
			TestSofeaLute::assertResults_Binding1CC8Mut3Flex4_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_MultiSweepLoMem() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut3Flex4_Standard.get(),
			50.0,
			16*1024,
			TestSofeaLute::assertResults_Binding1CC8Mut3Flex4_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_MultiSweepLoMem_2Threads() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut3Flex4_Standard.get(),
			50.0,
			16*1024,
			2,
			TestSofeaLute::assertResults_Binding1CC8Mut3Flex4_Standard
		);
	}
	@Test
	public void test_Binding1CC8Mut3Flex4_Standard_MultiSweepLoMem_4Threads() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut3Flex4_Standard.get(),
			50.0,
			16*1024,
			4,
			TestSofeaLute::assertResults_Binding1CC8Mut3Flex4_Standard
		);
	}


	@Test
	public void test_Binding1CC8Mut2Flex2_NoPLUG_LeafCounts() {
		Design design = Designs.Binding1CC8Mut2Flex2_NoPLUG.get();
		assertLeafCounts(design);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_NoPLUG_ZValueBounds() {
		Design design = Designs.Binding1CC8Mut2Flex2_NoPLUG.get();
		assertZValueBounds(design);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_NoPLUG_ZSumBounds() {
		Design design = Designs.Binding1CC8Mut2Flex2_NoPLUG.get();
		assertZSumBounds(design);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_NoPLUG_CalcG() {
		Design design = Designs.Binding1CC8Mut2Flex2_NoPLUG.get();
		assertGStates(design, Arrays.asList("SER", "GLY"), -17.356, -79.225, -49.195);
		assertGStates(design, Arrays.asList("VAL", "GLY"), -8.825, Double.POSITIVE_INFINITY, -49.195);
		assertGStates(design, Arrays.asList("SER", "VAL"), -12.050, Double.POSITIVE_INFINITY, -49.195);
		assertGStates(design, Arrays.asList("VAL", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, -49.195);
	}
	private static void assertResults_Binding1CC8Mut2Flex2_NoPLUG(Results results) {
		results.assertGUnsequenced(-49.195);
		results.assertGSequenced(Arrays.asList("SER", "GLY"), -17.356, -79.225);
		results.assertGSequenced(Arrays.asList("VAL", "GLY"), -8.825, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("SER", "VAL"), -12.050, Double.POSITIVE_INFINITY);
		results.assertGSequenced(Arrays.asList("VAL", "VAL"), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_NoPLUG_SingleSweep() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut2Flex2_NoPLUG.get(),
			Double.POSITIVE_INFINITY,
			1024*1024,
			TestSofeaLute::assertResults_Binding1CC8Mut2Flex2_NoPLUG
		);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_NoPLUG_MultiSweepHiMem() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut2Flex2_NoPLUG.get(),
			50.0,
			1024*1024,
			TestSofeaLute::assertResults_Binding1CC8Mut2Flex2_NoPLUG
		);
	}
	@Test
	public void test_Binding1CC8Mut2Flex2_NoPLUG_MultiSweepLoMem() {
		sweepUntilExhaustion(
			Designs.Binding1CC8Mut2Flex2_NoPLUG.get(),
			50.0,
			16*1024,
			TestSofeaLute::assertResults_Binding1CC8Mut2Flex2_NoPLUG
		);
	}


	/** brute forces every node in the tree and calls the supplied block with a ConfIndex instance describing the node */
	public static void forEachNode(SofeaLute.StateInfo stateInfo, Consumer<ConfIndex> block) {

		ConfIndex index = stateInfo.makeConfIndex();

		// NOTE: java has a hard time with recursive lambdas,
		// so use an array to work around the compiler's limitations
		Runnable[] f = { null };
		f[0] = () -> {

			// call the supplied block
			block.accept(index);

			// stop recursion if this is a leaf node
			if (index.isFullyDefined()) {
				return;
			}

			// otherwise, recurse
			int pos = index.numDefined;
			for (int rc : stateInfo.rcs.get(pos)) {
				index.assignInPlace(pos, rc);
				f[0].run();
				index.unassignInPlace(pos);
			}
		};
		f[0].run();
	}

	public static void assertLeafCounts(Design design) {

		SofeaLute sofea = new SofeaLute.Builder(design.confSpace)
			.configEachState(state -> design.config[state.index])
			.build();

		for (MultiStateConfSpace.State state : design.confSpace.states) {
			SofeaLute.StateInfo stateInfo = sofea.getStateInfo(state);

			forEachNode(stateInfo, index -> {
				BigIntegerBounds bounds = stateInfo.boundLeavesPerSequence(index);
				Map<Sequence,BigInteger> counts = stateInfo.countLeavesBySequence(index);
				BigInteger minCount = counts.values().stream().min(BigInteger::compareTo).orElse(null);
				BigInteger maxCount = counts.values().stream().max(BigInteger::compareTo).orElse(null);
				assertThat(bounds.lower, lessThanOrEqualTo(minCount));
				assertThat(bounds.upper, greaterThanOrEqualTo(maxCount));
			});
		}
	}

	public static void assertZValueBounds(Design design) {

		SofeaLute sofea = new SofeaLute.Builder(design.confSpace)
			.configEachState(state -> design.config[state.index])
			.build();

		RCTuple tripleTuple = new RCTuple(0, 0, 0, 0, 0, 0);

		for (MultiStateConfSpace.State state : design.confSpace.states) {
			SofeaLute.StateInfo stateInfo = sofea.getStateInfo(state);

			forEachNode(stateInfo, index -> {
				BigDecimalBounds bounds = new BigDecimalBounds(
					stateInfo.optimizeZ(index, tripleTuple, MathTools.Optimizer.Minimize),
					stateInfo.optimizeZ(index, tripleTuple, MathTools.Optimizer.Maximize)
				);
				BigDecimalBounds exact = stateInfo.exactBoundZ(index, tripleTuple);
				assertThat(bounds, isAbsoluteBound(exact, 1e-3));
			});
		}
	}

	public static void assertZSumBounds(Design design) {

		SofeaLute sofea = new SofeaLute.Builder(design.confSpace)
			.configEachState(state -> design.config[state.index])
			.build();

		RCTuple tripleTuple = new RCTuple(0, 0, 0, 0, 0, 0);

		for (MultiStateConfSpace.State state : design.confSpace.states) {
			SofeaLute.StateInfo stateInfo = sofea.getStateInfo(state);

			forEachNode(stateInfo, index -> {

				// skip leaf nodes
				if (index.isFullyDefined()) {
					return;
				}

				BigDecimalBounds bounds = stateInfo.boundZ(index, tripleTuple, BigDecimal.ONE);
				BigDecimal exact = stateInfo.calcZ(index, stateInfo.rcs, BigDecimal.ONE);
				if (bounds != null) {
					assertThat(bounds, isAbsoluteBound(exact, 1e-3));
				} else {
					assertThat(exact.doubleValue(), isAbsolutely(0.0, 1e-3));
				}
			});
		}
	}

	public static void main(String[] args) {

		// calc exact state g values for all designs
		for (Designs designId : Designs.values()) {

			Design design = designId.get();
			log("%s", designId);

			for (Sequence seq : design.confSpace.seqSpace.getSequences()) {
				log("Arrays.asList(%s), %s",
					Streams.joinToString(seq.assignments(), ", ", a -> "\"" + a.getResType().name + "\""),
					Streams.joinToString(
						calcGStatesAStar(design, seq),
						", ",
						(double d) -> String.format("%.3f", d)
					)
				);
			}
		}
	}

	public static double[] calcGStatesAStar(Design design, Sequence seq) {
		return design.confSpace.states.stream()
			.mapToDouble(state -> {

				SofeaLute.StateConfig config = design.config[state.index];
				RCs rcs = new RCs(
					seq.makeRCs(state.confSpace),
					config.pmat
				);
				ConfAStarTree astar = new ConfAStarTree.Builder(design.emats[state.index], rcs)
					.setLUTE(config.luteEcalc)
					.build();

				BigMath z = new BigMath(bcalc.mathContext).set(0.0);
				for (ConfSearch.ScoredConf conf : astar.nextConfs(Double.POSITIVE_INFINITY)) {
					z.add(bcalc.calcPrecise(conf.getScore()));
				}

				return bcalc.freeEnergyPrecise(z.get());
			})
			.toArray();
	}

	public void assertGStates(Design design, List<String> resTypes, double ... expectedG) {
		assertGStates(
			design,
			design.confSpace.seqSpace.makeSequence(resTypes),
			expectedG
		);
	}

	public void assertGStates(Design design, Sequence seq, double ... expectedG) {

		SofeaLute sofea = new SofeaLute.Builder(design.confSpace)
			.configEachState(state -> design.config[state.index])
			.build();

		double[] g = design.confSpace.states.stream()
			.mapToDouble(state -> sofea.calcG(state, seq))
			.toArray();

		assertThat(g, isAbsolutely(expectedG, epsilonG));
	}

	private static class Results {

		public final Design design;
		public final Map<Sequence,SeqDB.SeqInfo> sequenced = new HashMap<>();
		public final List<BigDecimalBounds> unsequenced = new ArrayList<>();

		public Results(Design design, SeqDB seqdb) {

			this.design = design;

			for (Sequence seq : design.confSpace.seqSpace.getSequences()) {
				sequenced.put(seq, seqdb.getSequencedZSumBounds(seq));
			}
			for (MultiStateConfSpace.State state : design.confSpace.unsequencedStates) {
				unsequenced.add(seqdb.getUnsequencedZSumBounds(state));
			}
		}

		public void assertGUnsequenced(double ... expectedG) {
			assertThat(design.confSpace.unsequencedStates.size(), is(expectedG.length));
			for (MultiStateConfSpace.State state : design.confSpace.unsequencedStates) {
				double expG = expectedG[state.unsequencedIndex];
				DoubleBounds obsG = bcalc.freeEnergyPrecise(unsequenced.get(state.unsequencedIndex));
				assertThat(obsG, isAbsoluteBound(expG, epsilonG));
			}
		}

		public void assertGSequenced(List<String> resTypes, double ... expectedG) {
			assertThat(design.confSpace.sequencedStates.size(), is(expectedG.length));
			Sequence seq = design.confSpace.seqSpace.makeSequence(resTypes);
			SeqDB.SeqInfo seqInfo = sequenced.get(seq);
			for (MultiStateConfSpace.State state : design.confSpace.sequencedStates) {
				double expG = expectedG[state.sequencedIndex];
				DoubleBounds obsG = bcalc.freeEnergyPrecise(seqInfo.get(state));
				assertThat(obsG, isAbsoluteBound(expG, epsilonG));
			}
		}
	}

	private interface IntermediateChecker {
		void check(Results results);
	}

	public void sweepUntilExhaustion(Design design, double sweepDivisor, long fringeDBBytes, IntermediateChecker checker) {
		sweepUntilExhaustion(design, sweepDivisor, fringeDBBytes, 1, checker);
	}

	public void sweepUntilExhaustion(Design design, double sweepDivisor, long fringeDBBytes, int numThreads, IntermediateChecker checker) {
		try (TempFile fringedbFile = new TempFile(tmpdir, "fringe.db")) {
			try (TempFile seqdbFile = new TempFile(tmpdir, "seq.db")) {

				SofeaLute sofea = new SofeaLute.Builder(design.confSpace)
					.setFringeDBFile(fringedbFile)
					.setFringeDBBytes(fringeDBBytes)
					.setSeqDBFile(seqdbFile)
					.setSweepDivisor(sweepDivisor)
					.configEachState(state -> design.config[state.index])
					.setParallelism(Parallelism.makeCpu(numThreads))
					.build();

				sofea.init();

				// refine, and check results between each sweep
				sofea.refine((seqdb, fringedb, sweepCount) -> {
					checker.check(new Results(design, seqdb));
					return SofeaLute.Criterion.KeepIterating;
				});

				// check results once more at end, just for good measure
				try (SeqDB seqdb = sofea.openSeqDB()) {
					checker.check(new Results(design, seqdb));
				}
			}
		}
	}

	private static enum ConfSpaces {

		// WARNING: if you change any details here, make sure to delete any cached files in `TestSofea.tmpdir`

		Binding1CC8Flex4 {

			@Override
			public MultiStateConfSpace make() {

				Molecule pdb = PDBIO.readResource("/1CC8.ss.pdb");

				Strand design = new Strand.Builder(pdb)
					.setResidues("A68", "A73")
					.build();
				for (String resNum : Arrays.asList("A69", "A70")) {
					design.flexibility.get(resNum)
						.setLibraryRotamers(Strand.WildType)
						.addWildTypeRotamers()
						.setContinuous();
				}

				Strand target = new Strand.Builder(pdb)
					.setResidues("A2", "A67")
					.build();
				for (String resNum : Arrays.asList("A5", "A6")) {
					target.flexibility.get(resNum)
						.setLibraryRotamers(Strand.WildType)
						.addWildTypeRotamers()
						.setContinuous();
				}

				// make a multi-state conf space
				return new MultiStateConfSpace
					.Builder("design", new SimpleConfSpace.Builder().addStrands(design).build())
					.addMutableState("complex", new SimpleConfSpace.Builder().addStrands(design, target).build())
					.addUnmutableState("target", new SimpleConfSpace.Builder().addStrands(target).build())
					.build();
			}
		},

		Stability1CC8Mut2Flex2 {

			@Override
			public MultiStateConfSpace make() {

				Molecule pdb = PDBIO.readResource("/1CC8.ss.pdb");

				Strand design = new Strand.Builder(pdb)
					.setResidues("A68", "A73")
					.build();
				for (String resNum : Arrays.asList("A69", "A70")) {
					design.flexibility.get(resNum)
						.setLibraryRotamers(Strand.WildType, "VAL")
						.addWildTypeRotamers()
						.setContinuous();
				}

				// make a multi-state conf space
				return new MultiStateConfSpace
					.Builder("design", new SimpleConfSpace.Builder().addStrands(design).build())
					.build();
			}
		},

		Binding1CC8Mut2Flex2 {

			@Override
			public MultiStateConfSpace make() {

				Molecule pdb = PDBIO.readResource("/1CC8.ss.pdb");

				Strand design = new Strand.Builder(pdb)
					.setResidues("A68", "A73")
					.build();
				for (String resNum : Arrays.asList("A69", "A70")) {
					design.flexibility.get(resNum)
						.setLibraryRotamers(Strand.WildType, "VAL")
						.addWildTypeRotamers()
						.setContinuous();
				}

				Strand target = new Strand.Builder(pdb)
					.setResidues("A2", "A67")
					.build();
				for (String resNum : Arrays.asList("A5", "A6")) {
					target.flexibility.get(resNum)
						.setLibraryRotamers(Strand.WildType)
						.addWildTypeRotamers()
						.setContinuous();
				}

				// make a multi-state conf space
				return new MultiStateConfSpace
					.Builder("design", new SimpleConfSpace.Builder().addStrands(design).build())
					.addMutableState("complex", new SimpleConfSpace.Builder().addStrands(design, target).build())
					.addUnmutableState("target", new SimpleConfSpace.Builder().addStrands(target).build())
					.build();
			}
		},

		Binding1CC8Mut3Flex4 {

			@Override
			public MultiStateConfSpace make() {

				Molecule pdb = PDBIO.readResource("/1CC8.ss.pdb");

				Strand design = new Strand.Builder(pdb)
					.setResidues("A68", "A73")
					.build();
				for (String resNum : Arrays.asList("A69", "A70", "A71")) {
					design.flexibility.get(resNum)
						.setLibraryRotamers(Strand.WildType, "VAL")
						.addWildTypeRotamers()
						.setContinuous();
				}

				Strand target = new Strand.Builder(pdb)
					.setResidues("A2", "A67")
					.build();
				for (String resNum : Arrays.asList("A5", "A6", "A7", "A8")) {
					target.flexibility.get(resNum)
						.setLibraryRotamers(Strand.WildType)
						.addWildTypeRotamers()
						.setContinuous();
				}

				// make a multi-state conf space
				return new MultiStateConfSpace
					.Builder("design", new SimpleConfSpace.Builder().addStrands(design).build())
					.addMutableState("complex", new SimpleConfSpace.Builder().addStrands(design, target).build())
					.addUnmutableState("target", new SimpleConfSpace.Builder().addStrands(target).build())
					.build();
			}
		};

		private static final Map<ConfSpaces,MultiStateConfSpace> confSpaces = new EnumMap<>(ConfSpaces.class);

		public abstract MultiStateConfSpace make();

		public MultiStateConfSpace get() {
			return confSpaces.computeIfAbsent(this, key -> make());
		}
	}

	private static enum Designs {

		// WARNING: if you change any details here, make sure to delete any cached files in `TestSofea.tmpdir`

		Binding1CC8Flex4_Standard {
			@Override
			public Design make() {
				return Design.makeStandard(this, ConfSpaces.Binding1CC8Flex4.get());
			}
		},
		Stability1CC8Mut2Flex2_Standard {
			@Override
			public Design make() {
				return Design.makeStandard(this, ConfSpaces.Stability1CC8Mut2Flex2.get());
			}
		},
		Binding1CC8Mut2Flex2_Standard {
			@Override
			public Design make() {
				return Design.makeStandard(this, ConfSpaces.Binding1CC8Mut2Flex2.get());
			}
		},
		Binding1CC8Mut3Flex4_Standard {
			@Override
			public Design make() {
				return Design.makeStandard(this, ConfSpaces.Binding1CC8Mut3Flex4.get());
			}
		},
		Binding1CC8Mut2Flex2_NoPLUG {
			@Override
			public Design make() {
				return Design.makeNoPLUG(this, ConfSpaces.Binding1CC8Mut2Flex2.get());
			}
		};

		private static final Map<Designs,Design> designs = new EnumMap<>(Designs.class);

		public abstract Design make();

		public Design get() {
			return designs.computeIfAbsent(this, key -> make());
		}
	}

	private static class Design {

		public final MultiStateConfSpace confSpace;
		public final EnergyMatrix[] emats;
		public final SofeaLute.StateConfig[] config;

		public Design(MultiStateConfSpace confSpace) {
			this.confSpace = confSpace;
			this.emats = new EnergyMatrix[confSpace.states.size()];
			this.config = new SofeaLute.StateConfig[confSpace.states.size()];
		}

		private static interface StateConfigurator {
			void config(MultiStateConfSpace.State state, ConfEnergyCalculator confEcalc, Design design);
		}

		public static Design makeStandard(Designs id, MultiStateConfSpace confSpace) {
			return make(id, confSpace, (state, confEcalc, design) -> {

				EnergyMatrix emat = calcEmat(id, state, confEcalc);
				PruningMatrix pmat = calcPmat(id, state, emat, 50.0, true);
				LUTEState luteState = trainLute(id, state, confEcalc, emat, pmat);

				design.emats[state.index] = emat;
				design.config[state.index] = new SofeaLute.StateConfig(
					new LUTEConfEnergyCalculator(state.confSpace, luteState),
					pmat
				);
			});
		}

		public static Design makeNoPLUG(Designs id, MultiStateConfSpace confSpace) {
			return make(id, confSpace, (state, confEcalc, design) -> {

				EnergyMatrix emat = calcEmat(id, state, confEcalc);
				PruningMatrix pmat = calcPmat(id, state, emat, 10.0, false);
				LUTEState luteState = trainLute(id, state, confEcalc, emat, pmat);

				design.emats[state.index] = emat;
				design.config[state.index] = new SofeaLute.StateConfig(
					new LUTEConfEnergyCalculator(state.confSpace, luteState),
					pmat
				);
			});
		}

		private static Design make(Designs id, MultiStateConfSpace confSpace, StateConfigurator configurator) {

			Design design = new Design(confSpace);

			try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
				.setParallelism(fullCPUParallelism)
				.build()) {

				for (MultiStateConfSpace.State state : confSpace.states) {
					ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc).build();
					configurator.config(state, confEcalc, design);
				}
			}

			return design;
		}

		private static EnergyMatrix calcEmat(Designs id, MultiStateConfSpace.State state, ConfEnergyCalculator confEcalc) {
			return new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File(tmpdir, String.format("%s.%s.emat", id, state.name)))
				.build()
				.calcEnergyMatrix();
		}

		private static PruningMatrix calcPmat(Designs id, MultiStateConfSpace.State state, EnergyMatrix emat, double goldsteinThreshold, boolean doPLUG) {

			SimpleDEE.Runner runner = new SimpleDEE.Runner()
				.setCacheFile(new File(tmpdir, String.format("%s.%s.pmat", id, state.name)))
				.setParallelism(fullCPUParallelism)
				.setTransitivePruning(true) // NOTE: transitive pruning very important for LUTE
				.setShowProgress(true);

			runner.setThreshold(100.0);
			runner.setGoldsteinDiffThreshold(goldsteinThreshold);

			if (doPLUG) {
				runner.setSinglesPlugThreshold(0.6);
				runner.setPairsPlugThreshold(0.6);
				//runner.setTriplesPlugThreshold(0.6)  // PLUG triples take too long!
			}

			return runner.run(state.confSpace, emat);
		}

		private static LUTEState trainLute(Designs id, MultiStateConfSpace.State state, ConfEnergyCalculator confEcalc, EnergyMatrix emat, PruningMatrix pmat) {

			LUTEState luteState;

			File luteFile = new File(tmpdir, String.format("%s.%s.lute", id, state.name));
			if (luteFile.exists()) {

				luteState = LUTEIO.read(luteFile);
				log("read LUTE state from file: %s", luteFile.getAbsolutePath());

			} else {

				try (ConfDB confdb = new ConfDB(state.confSpace)) {
					ConfDB.ConfTable confTable = confdb.new ConfTable("lute");

					final int randomSeed = 12345;
					final LUTE.Fitter fitter = LUTE.Fitter.OLSCG;
					final double maxOverfittingScore = 1.5;
					final double maxRMSE = 0.1;

					// compute LUTE fit
					LUTE lute = new LUTE(state.confSpace);
					ConfSampler sampler = new RandomizedDFSConfSampler(state.confSpace, pmat, randomSeed);
					lute.sampleTuplesAndFit(confEcalc, emat, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE);
					lute.reportConfSpaceSize(pmat);

					luteState = new LUTEState(lute.getTrainingSystem());
					LUTEIO.write(luteState, luteFile);
					log("wrote LUTE state to file: %s", luteFile.getAbsolutePath());
				}
			}

			return luteState;
		}
	}
}
