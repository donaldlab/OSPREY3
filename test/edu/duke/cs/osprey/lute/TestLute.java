/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.lute;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.kstar.*;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Streams;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static edu.duke.cs.osprey.tools.Log.log;


public class TestLute {
	
	private static SimpleConfSpace protein;
	private static SimpleConfSpace ligand;
	private static SimpleConfSpace complex;
	private static ForcefieldParams ffparams;
	
	@BeforeClass
	public static void beforeClass() {

		// configure the forcefield
		ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		TestLute.protein = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		TestLute.ligand = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();
		TestLute.complex = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.build();
	}

	private static ConfEnergyCalculator makeConfEcalc(SimpleConfSpace confSpace, EnergyCalculator ecalc) {
		return new ConfEnergyCalculator.Builder(confSpace, ecalc)
			.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcReferenceEnergies()
			).build();
	}

	private static EnergyMatrix calcEmat(ConfEnergyCalculator confEcalc) {
		return new SimplerEnergyMatrixCalculator.Builder(confEcalc)
			.build()
			.calcEnergyMatrix();
	}

	private static PruningMatrix calcPmat(SimpleConfSpace confSpace, EnergyMatrix emat) {
		return new SimpleDEE.Runner()
			.setSinglesThreshold(100.0)
			.setPairsThreshold(100.0)
			.setGoldsteinDiffThreshold(100.0)
			.setShowProgress(true)
			.run(confSpace, emat);
	}
	
	public static void main(String[] args) {

		// do all the designs with regular machinery, to compute the "correct" answers

		beforeClass();

		// do an A* enumeration
		{
			SimpleConfSpace confSpace = complex;

			try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
				.setParallelism(Parallelism.makeCpu(4))
				.build()) {

				ConfEnergyCalculator confEcalc = makeConfEcalc(confSpace, ecalc);
				EnergyMatrix emat = calcEmat(confEcalc);

				// do a GMEC search
				ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
					.setTraditional()
					.build();

				List<ConfSearch.ScoredConf> confs = astar.nextConfs(-20.0);
				List<ConfSearch.EnergiedConf> econfs = confEcalc.calcAllEnergies(confs);
				econfs.sort(Comparator.comparing((conf) -> conf.getEnergy()));

				for (ConfSearch.EnergiedConf conf : econfs) {
					log("assertConf(astar.nextConf(), new int[] { %s }, %.6f, epsilon);",
						String.join(",", IntStream.of(conf.getAssignments())
							.mapToObj(i -> Integer.toString(i))
							.collect(Collectors.toList())
						),
						conf.getEnergy()
					);
				}
			}
		}

		// do a GMEC design on just the protein strand
		{
			SimpleConfSpace confSpace = protein;

			try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build()) {

				ConfEnergyCalculator confEcalc = makeConfEcalc(confSpace, ecalc);
				EnergyMatrix emat = calcEmat(confEcalc);

				// do a GMEC search
				ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
					.setTraditional()
					.build();
				Queue.FIFO<ConfSearch.EnergiedConf> confs = new SimpleGMECFinder.Builder(astar, confEcalc)
					.build()
					.find(1.0);

				log("GMEC confs:");
				log("assertThat(confs.size(), is(%dL));", confs.size());
				while (!confs.isEmpty()) {
					ConfSearch.EnergiedConf conf = confs.poll();
					log("assertConf(confs.poll(), new int[] { %s }, %.6f, epsilon);",
						String.join(",", IntStream.of(conf.getAssignments())
							.mapToObj(i -> Integer.toString(i))
							.collect(Collectors.toList())
						),
						conf.getEnergy()
					);
				}
			}
		}

		final double epsilon = 0.01;

		// compute the pfuncs of the wild-type sequence
		{
			try (EnergyCalculator ecalc = new EnergyCalculator.Builder(complex, ffparams).build()) {

				Function<SimpleConfSpace,PartitionFunction.Result> pfuncCalculator = (confSpace) -> {

					ConfEnergyCalculator confEcalc = makeConfEcalc(confSpace, ecalc);
					EnergyMatrix emat = calcEmat(confEcalc);

					// pick the wild-type sequence
					Sequence sequence = confSpace.makeWildTypeSequence();
					RCs rcs = sequence.makeRCs(confSpace);
					ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build();

					// estimate the pfunc
					GradientDescentPfunc pfunc = new GradientDescentPfunc(confEcalc);
					pfunc.init(astar, rcs.getNumConformations(), epsilon);
					pfunc.setStabilityThreshold(null);
					pfunc.compute();

					return pfunc.makeResult();
				};

				PartitionFunction.Result resultProtein = pfuncCalculator.apply(protein);
				log("protein wild-type pfunc: %s %.4f", resultProtein.toString(), resultProtein.values.getEffectiveEpsilon());
				PartitionFunction.Result resultLigand = pfuncCalculator.apply(ligand);
				log("ligand  wild-type pfunc: %s %.4f", resultLigand.toString(), resultLigand.values.getEffectiveEpsilon());
				PartitionFunction.Result resultComplex = pfuncCalculator.apply(complex);
				log("complex wild-type pfunc: %s %.4f", resultComplex.toString(), resultComplex.values.getEffectiveEpsilon());
			}
		}

		KStarScoreWriter.Formatter testFormatter = (KStarScoreWriter.ScoreInfo info) ->
			String.format("assertKstar(scoredSequences.get(%d), \"%s\", %.6f, %.6f, fudge); // protein %s ligand %s complex %s",
				info.sequenceNumber,
				info.sequence.toString(Sequence.Renderer.ResTypeMutations),
				info.kstarScore.lowerBoundLog10(),
				info.kstarScore.upperBoundLog10(),
				info.kstarScore.protein.toString(),
				info.kstarScore.ligand.toString(),
				info.kstarScore.complex.toString()
			);

		// do a K* design
		{
			try (EnergyCalculator ecalc = new EnergyCalculator.Builder(complex, ffparams).build()) {

				KStar.Settings settings = new KStar.Settings.Builder()
					.setEpsilon(epsilon)
					.setStabilityThreshold(null)
					.addScoreConsoleWriter(testFormatter)
					.build();
				KStar kstar = new KStar(protein, ligand, complex, settings);
				for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

					info.confEcalc = makeConfEcalc(info.confSpace, ecalc);
					EnergyMatrix emat = calcEmat(info.confEcalc);
					info.confSearchFactory = (rcs) ->
						new ConfAStarTree.Builder(emat, rcs)
							.setTraditional()
							.build();
				}
				kstar.run();
			}
		}

		// do a BBK* design
		{
			try (EnergyCalculator ecalc = new EnergyCalculator.Builder(complex, ffparams).build()) {

				// configure BBK*
				KStar.Settings kstarSettings = new KStar.Settings.Builder()
					.setEpsilon(epsilon)
					.setStabilityThreshold(null)
					.addScoreConsoleWriter(testFormatter)
					.build();
				BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
					.setNumBestSequences(5)
					.build();
				BBKStar bbkstar = new BBKStar(protein, ligand, complex, kstarSettings, bbkstarSettings);
				for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

					// minimized energies
					info.confEcalcMinimized = makeConfEcalc(info.confSpace, ecalc);
					EnergyMatrix emat = calcEmat(info.confEcalcMinimized);
					info.confSearchFactoryMinimized = (rcs) ->
						new ConfAStarTree.Builder(emat, rcs)
							.setTraditional()
							.build();

					// rigid energies
					EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalc)
						.setIsMinimizing(false)
						.build();
					ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(info.confEcalcMinimized, ecalcRigid);
					EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
						.build()
						.calcEnergyMatrix();
					info.confSearchFactoryRigid = (rcs) ->
						new ConfAStarTree.Builder(ematRigid, rcs)
							.setTraditional()
							.build();
				}
				bbkstar.run();
			}
		}
	}

	private static void assertConf(ConfSearch.ScoredConf conf, int[] assignments, double score, double epsilon) {
		assertThat(conf.getAssignments(), is(assignments));
		assertThat(conf.getScore(), isAbsolutely(score, epsilon));
	}

	private static LUTEConfEnergyCalculator train(SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc, EnergyMatrix emat, PruningMatrix pmat) {

		try (ConfDB confdb = new ConfDB(confSpace)) {
			ConfDB.ConfTable confTable = confdb.new ConfTable("lute");

			final int randomSeed = 12345;
			final LUTE.Fitter fitter = LUTE.Fitter.OLSCG;
			final double maxOverfittingScore = 1.5;
			final double maxRMSE = 0.1;

			// compute LUTE fit
			LUTE lute = new LUTE(confSpace);
			ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, randomSeed);
			lute.sampleTuplesAndFit(confEcalc, emat, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE);
			lute.reportConfSpaceSize(pmat);

			return new LUTEConfEnergyCalculator(confSpace, new LUTEState(lute.getTrainingSystem()));
		}
	}

	@Test
	public void astar() {

		SimpleConfSpace confSpace = complex;

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			ConfEnergyCalculator confEcalc = makeConfEcalc(confSpace, ecalc);
			EnergyMatrix emat = calcEmat(confEcalc);
			PruningMatrix pmat = calcPmat(confSpace, emat);

			// train LUTE
			LUTEConfEnergyCalculator luteEcalc = train(confSpace, confEcalc, emat, pmat);

			ConfAStarTree astar = new ConfAStarTree.Builder(null, pmat)
				.setLUTE(luteEcalc)
				.build();

			final double epsilon = 1e-1;
			assertConf(astar.nextConf(), new int[] { 13,13,11 }, -29.037828, epsilon);
			assertConf(astar.nextConf(), new int[] { 13,13,40 }, -28.836836, epsilon);
			assertConf(astar.nextConf(), new int[] { 11,13,11 }, -28.321740, epsilon);
			assertConf(astar.nextConf(), new int[] { 11,13,40 }, -28.152335, epsilon);
			assertConf(astar.nextConf(), new int[] { 13,13,9 }, -27.791179, epsilon);
			assertConf(astar.nextConf(), new int[] { 13,13,10 }, -27.072290, epsilon);
		}
	}

	@Test
	public void gmec() {

		// do a GMEC design on just one strand
		SimpleConfSpace confSpace = protein;

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			ConfEnergyCalculator confEcalc = makeConfEcalc(confSpace, ecalc);
			EnergyMatrix emat = calcEmat(confEcalc);
			PruningMatrix pmat = calcPmat(confSpace, emat);

			// train LUTE
			LUTEConfEnergyCalculator luteEcalc = train(confSpace, confEcalc, emat, pmat);

			// do a GMEC search
			Queue.FIFO<ConfSearch.ScoredConf> confs = new LUTEGMECFinder.Builder(pmat, luteEcalc)
				.build()
				.find(1.0);

			// LUTE should match two-position designs exactly
			assertThat(confs.size(), is(5L));
			final double epsilon = 1e-6;
			assertConf(confs.poll(), new int[] { 11,13 }, -3.713187, epsilon);
			assertConf(confs.poll(), new int[] { 12,13 }, -3.393738, epsilon);
			assertConf(confs.poll(), new int[] { 8,13 }, -3.170271, epsilon);
			assertConf(confs.poll(), new int[] { 13,13 }, -2.857468, epsilon);
			assertConf(confs.poll(), new int[] { 3,13 }, -2.737019, epsilon);
		}
	}

	@Test
	public void pfunc() {

		final double epsilon = 0.01;

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(complex, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			Function<SimpleConfSpace,PartitionFunction.Result> pfuncCalculator = (confSpace) -> {

				ConfEnergyCalculator confEcalc = makeConfEcalc(confSpace, ecalc);
				EnergyMatrix emat = calcEmat(confEcalc);
				PruningMatrix pmat = calcPmat(confSpace, emat);

				// train LUTE
				LUTEConfEnergyCalculator luteEcalc = train(confSpace, confEcalc, emat, pmat);

				// pick the wild-type sequence
				Sequence sequence = confSpace.makeWildTypeSequence();
				RCs rcs = new RCs(sequence.makeRCs(confSpace), pmat);
				ConfAStarTree astar = new ConfAStarTree.Builder(null, rcs)
					.setLUTE(luteEcalc)
					.build();

				// estimate the pfunc
				LUTEPfunc pfunc = new LUTEPfunc(luteEcalc);
				pfunc.init(astar, rcs.getNumConformations(), epsilon);
				pfunc.setStabilityThreshold(null);
				pfunc.setReportProgress(true);
				pfunc.compute();

				return pfunc.makeResult();
			};

			// LUTE energies should be exact here
			assertPfunc(pfuncCalculator.apply(protein), epsilon, 2.403883, 2.404058, 0.0);
			assertPfunc(pfuncCalculator.apply(ligand), epsilon, 11.209861, 11.209861, 0.0);
			// LUTE approximates energies here, so it doesn't compute quite the same pfunc values
			final double fudge = 0.01;
			assertPfunc(pfuncCalculator.apply(complex), epsilon, 21.553341, 21.554117, fudge);
		}
	}

	private static void assertPfunc(PartitionFunction.Result result, double epsilon, double lower, double upper, double fudge) {
		assertThat(result.values.getEffectiveEpsilon(), lessThanOrEqualTo(epsilon));
		assertThat(KStarScore.scoreToLog10(result.values.calcLowerBound()), lessThanOrEqualTo(lower*(1+fudge)));
		assertThat(KStarScore.scoreToLog10(result.values.calcUpperBound()), greaterThanOrEqualTo(upper*(1-fudge)));
	}

	@Test
	public void kstar() {

		final double epsilon = 0.01;

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(complex, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setStabilityThreshold(null)
				.build();
			KStar kstar = new KStar(protein, ligand, complex, settings);
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				info.confEcalc = makeConfEcalc(info.confSpace, ecalc);
				EnergyMatrix emat = calcEmat(info.confEcalc);
				PruningMatrix pmat = calcPmat(info.confSpace, emat);

				// train LUTE
				LUTEConfEnergyCalculator luteEcalc = train(info.confSpace, info.confEcalc, emat, pmat);

				info.confSearchFactory = (rcs) ->
					new ConfAStarTree.Builder(null, new RCs(rcs, pmat))
						.setLUTE(luteEcalc)
						.build();
			}

			List<KStar.ScoredSequence> scoredSequences = kstar.run();

			// LUTE approximates energies, so it doesn't compute quite the same K* scores
			final double fudge = 0.01;
			assertKstar(scoredSequences.get(0), "asp glu lys", 7.937791, 7.947125, fudge); // protein [2.402387 , 2.404539] (log10)                    ligand [11.206188,11.209861] (log10)                    complex [21.552191,21.555700] (log10)
			assertKstar(scoredSequences.get(1), "asp glu ASP", 3.438280, 3.444701, fudge); // protein [2.402387 , 2.404539] (log10)                    ligand [1.189648 , 1.190554] (log10)                    complex [7.033372 , 7.036735] (log10)
			assertKstar(scoredSequences.get(2), "asp glu GLU", 2.601122, 2.605583, fudge); // protein [2.402387 , 2.404539] (log10)                    ligand [1.247908 , 1.247908] (log10)                    complex [6.253568 , 6.255878] (log10)
			assertKstar(scoredSequences.get(3), "asp ASP lys", 6.068227, 6.077714, fudge); // protein [-1.111587,-1.108619] (log10)                    ligand [11.206188,11.209861] (log10)                    complex [16.169469,16.172315] (log10)
			assertKstar(scoredSequences.get(4), "GLU glu lys", 6.762618, 6.774455, fudge); // protein [3.069510 , 3.073791] (log10)                    ligand [11.206188,11.209861] (log10)                    complex [21.046269,21.050153] (log10)
		}
	}

	private static void assertKstar(KStar.ScoredSequence scoredSequence, String expectedSequence, double kstarLower, double kstarUpper, double fudge) {
		assertThat(scoredSequence.sequence.toString(Sequence.Renderer.ResTypeMutations), is(expectedSequence));
		assertThat(scoredSequence.score.lowerBoundLog10(), lessThanOrEqualTo(kstarLower*(1+fudge)));
		assertThat(scoredSequence.score.upperBoundLog10(), greaterThanOrEqualTo(kstarUpper*(1-fudge)));
	}

	@Test
	public void bbkstar() {

		final double epsilon = 0.01;

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(complex, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			KStar.Settings kstarSettings = new KStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setStabilityThreshold(null)
				.addScoreConsoleWriter()
				.build();
			BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
				.setNumBestSequences(5)
				.build();
			BBKStar bbkstar = new BBKStar(protein, ligand, complex, kstarSettings, bbkstarSettings);
			for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

				info.confEcalcMinimized = makeConfEcalc(info.confSpace, ecalc);
				EnergyMatrix emat = calcEmat(info.confEcalcMinimized);
				PruningMatrix pmat = calcPmat(info.confSpace, emat);

				// train LUTE
				LUTEConfEnergyCalculator luteEcalc = train(info.confSpace, info.confEcalcMinimized, emat, pmat);

				info.confSearchFactoryMinimized = (rcs) ->
					new ConfAStarTree.Builder(null, new RCs(rcs, pmat))
						.setLUTE(luteEcalc)
						.build();
				info.confSearchFactoryRigid = info.confSearchFactoryMinimized;
				info.pfuncFactory = new PartitionFunctionFactory(info.confSpace, info.confEcalcMinimized, info.id);
			}

			List<KStar.ScoredSequence> scoredSequences = bbkstar.run();

			// LUTE approximates energies, so it doesn't compute quite the same K* scores
			final double fudge = 0.01;
			assertKstar(scoredSequences.get(0), "asp glu lys", 7.937791, 7.947125, fudge); // protein [2.402387 , 2.404539] (log10)                    ligand [11.206188,11.209861] (log10)                    complex [21.552191,21.555700] (log10)
			assertKstar(scoredSequences.get(1), "GLU glu lys", 6.762618, 6.774455, fudge); // protein [3.069510 , 3.073791] (log10)                    ligand [11.206188,11.209861] (log10)                    complex [21.046269,21.050153] (log10)
			assertKstar(scoredSequences.get(2), "asp ASP lys", 6.068227, 6.077714, fudge); // protein [-1.111587,-1.108619] (log10)                    ligand [11.206188,11.209861] (log10)                    complex [16.169469,16.172315] (log10)
			assertKstar(scoredSequences.get(3), "asp glu ASP", 3.438280, 3.444701, fudge); // protein [2.402387 , 2.404539] (log10)                    ligand [1.189648 , 1.190554] (log10)                    complex [7.033372 , 7.036735] (log10)
			assertKstar(scoredSequences.get(4), "asp glu GLU", 2.601122, 2.605583, fudge); // protein [2.402387 , 2.404539] (log10)                    ligand [1.247908 , 1.247908] (log10)                    complex [6.253568 , 6.255878] (log10)
		}
	}

	@Test
	public void sampleCaching() {

		// every time we increase the size of the LUTE model (by adding tuples),
		// the sampler should reuse samples from the previous model size

		Strand protein = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20", "A21")) {
			protein.flexibility.get(resNum).setLibraryRotamers("VAL");
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrands(protein).build();

		// calc the emat
		EnergyMatrix emat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()
		) {
			emat = calcEmat(makeConfEcalc(confSpace, ecalc));
		}

		// use an empty pmat
		PruningMatrix pmat = new PruningMatrix(confSpace);

		final int randomSeed = 12345;
		ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, randomSeed);

		LUTE lute = new LUTE(confSpace);

		Supplier<Set<int[]>> getSamples = () -> {
			for (int i=0; i<3; i++) {
				sampler.sampleConfsForTuples(lute.trainingSet, i);
				sampler.sampleConfsForTuples(lute.testSet, i);
			}
			Conf.Set out = new Conf.Set();
			out.addAll(lute.trainingSet.getAllConfs());
			out.addAll(lute.testSet.getAllConfs());
			return out;
		};

		// samples confs for pair tuples
		lute.addTuples(lute.getUnprunedPairTuples(pmat));
		Set<int[]> pairsSamples = getSamples.get();

		// sample one round of sparse triples
		lute.addTuples(lute.sampleTripleTuplesByStrongInteractions(emat, pmat, 1));
		Set<int[]> triples1Samples = getSamples.get();

		// we should re-use all of the samples from the pairs
		assertThat(triples1Samples.containsAll(pairsSamples), is(true));

		// but the triples should have extra samples too
		assertThat(triples1Samples.size(), greaterThan(pairsSamples.size()));

		// sample another round of sparse triples
		lute.addUniqueTuples(lute.sampleTripleTuplesByStrongInteractions(emat, pmat, 2));
		Set<int[]> triples2Samples = getSamples.get();

		// we should re-use all of the samples from the previous sets
		assertThat(triples2Samples.containsAll(pairsSamples), is(true));
		assertThat(triples2Samples.containsAll(triples1Samples), is(true));

		// but the this one should have extra samples too
		assertThat(triples2Samples.size(), greaterThan(pairsSamples.size()));
		assertThat(triples2Samples.size(), greaterThan(triples1Samples.size()));
	}

	@Test
	public void sampleCachingAcrossRuns() {

		// if we're re-doing the same design with the same random seed,
		// make sure we sample confs in the same order

		Strand protein = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20", "A21")) {
			protein.flexibility.get(resNum).setLibraryRotamers("VAL");
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrands(protein).build();

		// calc the emat
		EnergyMatrix emat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()
		) {
			emat = calcEmat(makeConfEcalc(confSpace, ecalc));
		}

		// just use a blank pmat
		PruningMatrix pmat = new PruningMatrix(confSpace);

		// use the same seed for all LUTE runs
		int randomSeed = 12345;

		BiFunction<LUTE,ConfSampler,Set<int[]>> getSamples = (lute, sampler) -> {
			for (int i=0; i<3; i++) {
				sampler.sampleConfsForTuples(lute.trainingSet, i);
				sampler.sampleConfsForTuples(lute.testSet, i);
			}
			Conf.Set out = new Conf.Set();
			out.addAll(lute.trainingSet.getAllConfs());
			out.addAll(lute.testSet.getAllConfs());
			return out;
		};

		Function<LUTE,List<RCTuple>> getTuples = (lute) ->
			Streams.of(lute.tuplesIndex.iterator()).collect(Collectors.toList());

		// save the samples from the initial LUTE run
		final List<RCTuple> pairsTuples;
		final List<RCTuple> triples1Tuples;
		final List<RCTuple> triples2Tuples;
		final Set<int[]> pairsSamples;
		final Set<int[]> triples1Samples;
		final Set<int[]> triples2Samples;
		{
			LUTE lute = new LUTE(confSpace);
			ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, randomSeed);

			// start with just pair tuples
			lute.addUniqueTuples(lute.getUnprunedPairTuples(pmat));
			pairsTuples = getTuples.apply(lute);
			pairsSamples = getSamples.apply(lute, sampler);

			// then add triples
			lute.addUniqueTuples(lute.sampleTripleTuplesByStrongInteractions(emat, pmat, 1));
			triples1Tuples = getTuples.apply(lute);
			triples1Samples = getSamples.apply(lute, sampler);

			// then add more triples
			lute.addUniqueTuples(lute.sampleTripleTuplesByStrongInteractions(emat, pmat, 2));
			triples2Tuples = getTuples.apply(lute);
			triples2Samples = getSamples.apply(lute, sampler);
		}

		// do another LUTE run and check the samples
		{
			LUTE lute = new LUTE(confSpace);
			ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, randomSeed);

			// start with just pair tuples
			lute.addUniqueTuples(lute.getUnprunedPairTuples(pmat));
			assertThat(getTuples.apply(lute), is(pairsTuples));
			assertThat(getSamples.apply(lute, sampler), is(pairsSamples));

			// then add triples
			lute.addUniqueTuples(lute.sampleTripleTuplesByStrongInteractions(emat, pmat, 1));
			assertThat(getTuples.apply(lute), is(triples1Tuples));
			assertThat(getSamples.apply(lute, sampler), is(triples1Samples));

			// then add more triples
			lute.addUniqueTuples(lute.sampleTripleTuplesByStrongInteractions(emat, pmat, 2));
			assertThat(getTuples.apply(lute), is(triples2Tuples));
			assertThat(getSamples.apply(lute, sampler), is(triples2Samples));
		}
	}

	@Test
	public void forEachInAfterAddTriples() {

		// make sure adding tuples to the model works

		Strand protein = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19")) {
			protein.flexibility.get(resNum).setLibraryRotamers("VAL"); // 3 rotamers each
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrands(protein).build();

		// just use a blank pmat
		PruningMatrix pmat = new PruningMatrix(confSpace);

		// start with just pair tuples
		LUTE lute = new LUTE(confSpace);
		lute.addUniqueTuples(lute.getUnprunedPairTuples(pmat));

		// pick a conf, any conf
		int[] conf = new int[] { 1, 2, 0, 1 };

		// check the conf for pairs
		Set<RCTuple> tupleIndices = new HashSet<>();
		lute.tuplesIndex.forEachIn(conf, false, false, (t) ->
			tupleIndices.add(lute.tuplesIndex.get(t))
		);
		assertThat(tupleIndices, containsInAnyOrder(
			new RCTuple(0, 1, 1, 2).sorted(),
			new RCTuple(0, 1, 2, 0).sorted(),
			new RCTuple(0, 1, 3, 1).sorted(),
			new RCTuple(1, 2, 2, 0).sorted(),
			new RCTuple(1, 2, 3, 1).sorted(),
			new RCTuple(2, 0, 3, 1).sorted()
		));

		// add some triple tuples
		lute.addUniqueTuples(Arrays.asList(
			new RCTuple(0, 1, 1, 2, 2, 0).sorted(), // in the conf
			new RCTuple(1, 2, 2, 0, 3, 1).sorted(), // in the conf
			new RCTuple(0, 0, 1, 0, 2, 0).sorted() // not in the conf
		));

		// check the conf for pairs and triples
		tupleIndices.clear();
		lute.tuplesIndex.forEachIn(conf, false, false, (t) ->
			tupleIndices.add(lute.tuplesIndex.get(t))
		);
		assertThat(tupleIndices, containsInAnyOrder(
			new RCTuple(0, 1, 1, 2).sorted(),
			new RCTuple(0, 1, 2, 0).sorted(),
			new RCTuple(0, 1, 3, 1).sorted(),
			new RCTuple(1, 2, 2, 0).sorted(),
			new RCTuple(1, 2, 3, 1).sorted(),
			new RCTuple(2, 0, 3, 1).sorted(),
			new RCTuple(0, 1, 1, 2, 2, 0).sorted(),
			new RCTuple(1, 2, 2, 0, 3, 1).sorted()
		));
	}
}
