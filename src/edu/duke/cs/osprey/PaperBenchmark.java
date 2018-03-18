package edu.duke.cs.osprey;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.gpu.cuda.Gpus;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.BiFunction;

public class PaperBenchmark {

	public static final int NumTrials = 3;

	public static class Result {

		public class Trial {

			public final Stopwatch stopwatch = new Stopwatch();
			public double ops = 0.0;

			public Trial() {
				trials.add(this);
			}

			public Trial(Runnable task) {
				this();
				start();
				task.run();
				stop();
			}

			public Trial start() {
				stopwatch.start();
				return this;
			}

			public Trial stop() {
				stopwatch.stop();
				ops = (double)numRuns/stopwatch.getTimeS();
				return this;
			}
		}

		public final Parallelism parallelism;
		public final int numRuns;
		public final List<Trial> trials = new ArrayList<>();

		public Result(Parallelism parallelism, int numRuns) {
			this.parallelism = parallelism;
			this.numRuns = numRuns;
		}

		public double getBestOps() {
			return trials.stream()
				.map((trial) -> trial.ops)
				.max(Comparator.naturalOrder())
				.orElse(0.0);
		}
	}

	public static void main(String[] args)
		throws Exception {

		List<Parallelism> parallelisms = Arrays.asList(
			/*
			Parallelism.make(1, 0, 0),
			Parallelism.make(2, 0, 0),
			Parallelism.make(4, 0, 0),
			Parallelism.make(8, 0, 0),
			Parallelism.make(16, 0, 0),
			Parallelism.make(32, 0, 0),

			Parallelism.make(4, 1, 1),
			Parallelism.make(4, 1, 2),
			Parallelism.make(4, 1, 4),
			Parallelism.make(4, 1, 8),
			Parallelism.make(4, 1, 16),
			Parallelism.make(4, 1, 32)
			*/

			Parallelism.make(4, 4, 1), // 4
			Parallelism.make(4, 4, 2), // 8
			Parallelism.make(4, 4, 4), // 16
			Parallelism.make(4, 4, 8), // 32
			Parallelism.make(4, 4, 16), // 64
			Parallelism.make(4, 4, 32), // 128
			Parallelism.make(4, 4, 64) // 256
		);

		try (FileWriter out = new FileWriter(new File("benchmark.osprey3.tsv"))) {

			out.write("name\tnum CPU\tnum GPU\tum streams\tnum runs");
			for (int i=0; i<NumTrials; i++) {
				out.write("\tms " + i);
			}
			for (int i=0; i<NumTrials; i++) {
				out.write("\tops " + i);
			}
			out.write("\tbest ops\n");

			out.write("\n");
			//report("vsOsprey2", vsOsprey2(Parallelism.makeCpu(1)), out);

			out.write("\n");
			for (Parallelism parallelism : parallelisms) {
				report("pairMin", pairMin(parallelism), out);
			}

			out.write("\n");
			for (Parallelism parallelism : parallelisms) {
				report("confMin-1", confMin(1, 60, parallelism), out);
			}

			out.write("\n");
			for (Parallelism parallelism : parallelisms) {
				report("confMin-5", confMin(5, 10, parallelism), out);
			}

			out.write("\n");
			for (Parallelism parallelism : parallelisms) {
				report("confMin-10", confMin(10, 4, parallelism), out);
			}

			out.write("\n");
			for (Parallelism parallelism : parallelisms) {
				report("confMin-20", confMin(20, 2, parallelism), out);
			}
		}
	}

	public static void report(String name, Result result, FileWriter out)
		throws IOException {

		// report to the consle
		System.out.print(String.format("benchmark %30s %2d:%2d:%2d",
			name,
			result.parallelism.numThreads,
			result.parallelism.numGpus,
			result.parallelism.numStreamsPerGpu
		));
		if (result.numRuns == 0) {
			System.out.println("      skipped, lack of hardware");
		} else {
			for (Result.Trial trial : result.trials) {
				System.out.print(String.format("     %8s = %8.2f ops",
					trial.stopwatch.getTime(2),
					trial.ops
				));
			}
			System.out.println(String.format("     best ops: %.2f", result.getBestOps()));
		}

		// write to the file
		out.write(String.format("%s\t%d\t%d\t%d\t%d",
			name,
			result.parallelism.numThreads,
			result.parallelism.numGpus,
			result.parallelism.numStreamsPerGpu,
			result.numRuns
		));
		if (result.numRuns > 0) {
			for (Result.Trial trial : result.trials) {
				out.write(String.format("\t%.0f", trial.stopwatch.getTimeMs()));
			}
			for (Result.Trial trial : result.trials) {
				out.write(String.format("\t%.2f", trial.ops));
			}
			out.write(String.format("\t%.2f", result.getBestOps()));
		}
		out.write("\n");

		out.flush();
	}

	public static Result vsOsprey2(Parallelism parallelism) {

		Result result = new Result(parallelism, 1);
		Result.Trial trial = result.new Trial().start();

		Molecule mol = PDBIO.readFile("examples/paper-benchmark/1TP5_H_noHet_rename_final.pdb");

		Strand strand1 = new Strand.Builder(mol)
			.setResidues(302, 415)
			.build();
		strand1.flexibility.get(323).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand1.flexibility.get(325).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand1.flexibility.get(326).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		Strand strand2 = new Strand.Builder(mol)
			.setResidues(420, 425)
			.build();
		strand2.flexibility.get(423).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand2.flexibility.get(424).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand2.flexibility.get(425).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrands(strand1, strand2)
			.build();

		new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(parallelism)
			.use((ecalc) -> {

				SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcReferenceEnergies();

				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.setReferenceEnergies(eref)
					.build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
					.build()
					.calcEnergyMatrix();

				ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
					.setTraditional()
					.setShowProgress(true)
					.build();

				new SimpleGMECFinder.Builder(astar, confEcalc)
					.build()
					.find();
			});

		trial.stop();
		return result;
	}

	public static Result go1CC8(int numPos, Parallelism parallelism, BiFunction<SimpleConfSpace,ConfEnergyCalculator,Result> func) {

		// can we support the desired parallelism?
		int numCpuCores = Runtime.getRuntime().availableProcessors();
		if (parallelism.numThreads > numCpuCores) {
			return new Result(parallelism, 0);
		}
		if (parallelism.numGpus > Gpus.get().getGpus().size()) {
			return new Result(parallelism, 0);
		}

		Molecule mol = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");

		Strand strand = new Strand.Builder(mol).build();
		int firstResNum = 2;
		int lastResNum = firstResNum + numPos;
		for (int i=firstResNum; i<=lastResNum; i++) {
			strand.flexibility.get(i).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		AtomicReference<Result> resultRef = new AtomicReference<>();

		new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(parallelism)
			.use((ecalc) -> {

				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.build();

				resultRef.set(func.apply(confSpace, confEcalc));
			});

		return resultRef.get();
	}

	public static void calcEnergies(int numRuns, ConfEnergyCalculator confEcalc, RCTuple frag, ResidueInteractions inters) {
		for (int i=0; i<numRuns; i++) {
			confEcalc.calcEnergyAsync(frag, inters, (energy) -> {
				// ignore result, we're just benchmarking
			});
		}
		confEcalc.tasks.waitForFinish();
	}

	public static int calcNumAtomPairs(EnergyCalculator ecalc, SimpleConfSpace confSpace, RCTuple frag, ResidueInteractions inters) {
		ParametricMolecule mol = confSpace.makeMolecule(frag);
		ResidueForcefieldEnergy ff = new ResidueForcefieldEnergy(ecalc.resPairCache, inters, mol.mol);
		int numAtomPairs = 0;
		for (ResPairCache.ResPair resPair : ff.resPairs) {
			numAtomPairs += resPair.info.numAtomPairs;
		}
		return numAtomPairs;
	}

	public static Result warmupAndTime(int numWarmupsPerThread, int numRunsPerThread, ConfEnergyCalculator confEcalc, RCTuple frag, ResidueInteractions inters) {

		// warmup
		calcEnergies(numWarmupsPerThread*confEcalc.tasks.getParallelism(), confEcalc, frag, inters);

		// time it for reals
		Result result = new Result(confEcalc.ecalc.parallelism, numRunsPerThread*confEcalc.tasks.getParallelism());
		for (int i=0; i<NumTrials; i++) {
			result.new Trial(() -> calcEnergies(result.numRuns, confEcalc, frag, inters));
		}
		return result;
	}

	public static Result pairMin(Parallelism parallelism) {
		return go1CC8(2, parallelism, (confSpace, confEcalc) -> {

			RCTuple frag = new RCTuple(0, 0, 1, 0);
			ResidueInteractions inters = ResInterGen.of(confSpace)
				.addIntras(frag)
				.addInters(frag)
				.make();

			System.out.println("num atom pairs: " + calcNumAtomPairs(confEcalc.ecalc, confSpace, frag, inters));

			return warmupAndTime(100, 1000, confEcalc, frag, inters);
		});
	}

	public static Result confMin(int numPos, int numRunsPerThread, Parallelism parallelism) {
		return go1CC8(numPos, parallelism, (confSpace, confEcalc) -> {

			RCTuple frag = new RCTuple();
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				frag.pos.add(pos.index);
				frag.RCs.add(0);
			}
			ResidueInteractions inters = ResInterGen.of(confSpace)
				.addIntras(frag)
				.addInters(frag)
				.addShell(frag)
				.make();

			System.out.println("num atom pairs: " + calcNumAtomPairs(confEcalc.ecalc, confSpace, frag, inters));

			return warmupAndTime(1, numRunsPerThread, confEcalc, frag, inters);
		});
	}
}