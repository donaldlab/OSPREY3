package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


public class LUTELab {

	public static void main(String[] args) {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		List<String> mutableResNums = Arrays.asList("A38", "A39");
		//List<String> mutableResNums = Arrays.asList("A38", "A39", "A40"); // this works great on pairs
		//List<String> mutableResNums = Arrays.asList("A38", "A39", "A40", "A41"); // NOTE: need triples to get a good LUTE fit here
		for (String resNum : mutableResNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType, "VAL", "LEU")
				.addWildTypeRotamers()
				.setContinuous();
		}

		// add some flexibility
		List<String> flexibleResNums = Arrays.asList("A40", "A41", "A42", "A43", "A44");
		for (String resNum : flexibleResNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.make(8, 0, 0))
			//.setParallelism(Parallelism.make(8, 3, 16))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

			// compute conventional energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File("lute-test.emat.dat"))
				.build()
				.calcEnergyMatrix();

			// run DEE (super important for good LUTE fits!!)
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setSinglesThreshold(100.0)
				.setPairsThreshold(100.0)
				.setGoldsteinDiffThreshold(6.0)
				.setShowProgress(true)
				.run(confSpace, emat);
			//PruningMatrix pmat = new PruningMatrix(confSpace);

			Supplier<ConfSearch> astarFactory = () -> new ConfAStarTree.Builder(emat, pmat)
				.setTraditional()
				.build();

			LUTE lute;

			final File confDBFile = new File("lute-test.conf.db");
			try (ConfDB confdb = new ConfDB(confSpace, confDBFile)) {
				ConfDB.ConfTable confTable = confdb.new ConfTable("lute-test");

				// compute energies for some low-energy confs
				final int numConfs = 10;
				log("computing energies for %d confs...", numConfs);
				List<ConfSearch.ScoredConf> confs = new ArrayList<>();
				{
					ConfSearch astar = astarFactory.get();
					for (int i=0; i<numConfs; i++) {
						ConfSearch.ScoredConf conf = astar.nextConf();
						if (conf == null) {
							break;
						}
						confs.add(conf);
					}
				}
				List<ConfSearch.EnergiedConf> econfs = confEcalc.calcAllEnergies(confs, true, confTable);

				log("\nLUTE:\n");

				final int randomSeed = 12345;
				final int minSamplesPerTuple = 10;

				// compute LUTE matrix for pair tuples
				lute = new LUTE(confSpace);
				lute.addUnprunedPairTuples(pmat);
				lute.addUnprunedTripleTuples(pmat);
				// TODO: add sparse triples
				ConfSampler sampler = new UniformConfSampler(confSpace, randomSeed);
				/*ConfSampler sampler = new LowEnergyConfSampler(confSpace, randomSeed, pmat, (rcs) ->
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build()
				);*/
				confEcalc.resetCounters();
				lute.fit(confEcalc, confTable, sampler, 1.5);
				lute.reportConfSpaceSize(astarFactory.get());

				// compare conf energies
				LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(lute, ecalc);
				for (ConfSearch.EnergiedConf econf : econfs) {
					double luteEnergy = luteEcalc.calcEnergy(econf.getAssignments());
					log("conf %30s   score %9.4f      energy %9.4f   gap %7.4f      LUTE energy %9.4f   diff %7.4f",
						Arrays.toString(econf.getAssignments()),
						econf.getScore(),
						econf.getEnergy(),
						econf.getEnergy() - econf.getScore(),
						luteEnergy,
						luteEnergy - econf.getEnergy()
					);
				}

				final double pfuncEpsilon = 0.1;

				List<Sequence> sequences = getAllSequences(confSpace);
				log("sequences: %d", sequences.size());

				BiFunction<Sequence,PartitionFunction,PartitionFunction.Result> calcPfunc = (sequence, pfunc) -> {
					RCs rcs = sequence.makeRCs();
					ConfSearch astar = new ConfAStarTree.Builder(emat, new RCs(rcs, pmat))
						.setTraditional()
						.build();
					pfunc.init(astar, rcs.getNumConformations(), pfuncEpsilon);
					//pfunc.setReportProgress(true);
					pfunc.compute();
					return pfunc.makeResult();
				};

				// run traditional K*
				log("\nplain old K*:\n");
				confEcalc.resetCounters();
				Stopwatch oldSw = new Stopwatch().start();
				GradientDescentPfunc oldPfunc = new GradientDescentPfunc(confEcalc);
				oldPfunc.setConfTable(confTable);
				for (Sequence sequence : sequences) {
				//{ Sequence sequence = sequences.get(1);
					log("\t%s   %s   (%d energies)", sequence, calcPfunc.apply(sequence, oldPfunc), confEcalc.getNumRequests());
				}
				log("done in %s, %d energy calculations", oldSw.stop().getTime(2), confEcalc.getNumRequests());

				// use LUTE to do the same thing
				log("\nLUTE-powered K*:\n");
				Stopwatch luteSw = new Stopwatch().start();
				GradientDescentPfunc lutePfunc = new GradientDescentPfunc(luteEcalc);
				for (Sequence sequence : sequences) {
				//{ Sequence sequence = sequences.get(1);
					log("\t%s   %s   (%d energies)", sequence, calcPfunc.apply(sequence, lutePfunc), luteEcalc.getNumRequests());
				}
				log("done in %s, %d energy calculations", luteSw.stop().getTime(2), luteEcalc.getNumRequests());
			}
		}
	}

	private static List<Sequence> getAllSequences(SimpleConfSpace confSpace) {

		List<Sequence> sequences = new ArrayList<>();

		for (List<SimpleConfSpace.Position> mutablePositions : MathTools.powerset(confSpace.positions)) {

			// collect the mutations (res types except for wild type) for these positions into a simple list list
			List<List<String>> resTypes = new ArrayList<>();
			for (SimpleConfSpace.Position pos : mutablePositions) {
				resTypes.add(pos.resFlex.resTypes.stream()
					.filter((resType) -> !resType.equals(pos.resFlex.wildType))
					.collect(Collectors.toList())
				);
			}

			// enumerate all the combinations of res types
			for (List<String> mutations : MathTools.cartesianProduct(resTypes)) {

				// build the complex sequence
				Sequence complexSequence = confSpace.makeWildTypeSequence();
				for (int i=0; i<mutablePositions.size(); i++) {
					complexSequence.set(mutablePositions.get(i), mutations.get(i));
				}
				sequences.add(complexSequence);
			}
		}

		return sequences;
	}
}
