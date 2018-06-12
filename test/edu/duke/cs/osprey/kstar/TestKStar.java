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

package edu.duke.cs.osprey.kstar;

import static edu.duke.cs.osprey.TestBase.TempFile;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.Test;

import java.io.File;
import java.util.List;
import java.util.function.Function;


public class TestKStar {

	public static class ConfSpaces {
		public ForcefieldParams ffparams;
		public SimpleConfSpace protein;
		public SimpleConfSpace ligand;
		public SimpleConfSpace complex;
	}

	public static class Result {
		public KStar kstar;
		public List<KStar.ScoredSequence> scores;
	}

	public static Result runKStar(ConfSpaces confSpaces, double epsilon, String confDBPattern, boolean useExternalMemory, int maxSimultaneousMutations) {

		Parallelism parallelism = Parallelism.makeCpu(4);

		// how should we compute energies of molecules?
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
			.setParallelism(parallelism)
			.build()) {

			KStarScoreWriter.Formatter testFormatter = (KStarScoreWriter.ScoreInfo info) -> {

				Function<PartitionFunction.Result,String> formatPfunc = (pfuncResult) -> {
					if (pfuncResult.status == PartitionFunction.Status.Estimated) {
						return String.format("%12e", pfuncResult.values.qstar.doubleValue());
					}
					return "null";
				};

				return String.format("assertSequence(result, %3d, \"%s\", %-12s, %-12s, %-12s, epsilon); // protein %s ligand %s complex %s K* = %s",
					info.sequenceNumber,
					info.sequence.toString(Sequence.Renderer.ResType),
					formatPfunc.apply(info.kstarScore.protein),
					formatPfunc.apply(info.kstarScore.ligand),
					formatPfunc.apply(info.kstarScore.complex),
					info.kstarScore.protein.toString(),
					info.kstarScore.ligand.toString(),
					info.kstarScore.complex.toString(),
					info.kstarScore.toString()
				);
			};

			// configure K*
			KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setStabilityThreshold(null)
				.addScoreConsoleWriter(testFormatter)
				.setExternalMemory(useExternalMemory)
				.setMaxSimultaneousMutations(maxSimultaneousMutations)
				//.setShowPfuncProgress(true)
				.build();
			KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				// how should we define energies of conformations?
				info.confEcalc = new ConfEnergyCalculator.Builder(info.confSpace, ecalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(info.confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
					)
					.build();

				// calc energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.build()
					.calcEnergyMatrix();

				// how should confs be ordered and searched?
				info.confSearchFactory = (rcs) -> {
					ConfAStarTree.Builder builder = new ConfAStarTree.Builder(emat, rcs)
						.setTraditional();
					if (useExternalMemory) {
						builder.useExternalMemory();
					}
					return builder.build();
				};

				// set ConfDB if needed
				if (confDBPattern != null) {
					info.confDBFile = new File(confDBPattern.replace("*", info.type.name().toLowerCase()));
				}
			}

			// run K*
			Result result = new Result();
			result.kstar = kstar;
			result.scores = kstar.run();
			return result;
		}
	}

	public static ConfSpaces make2RL0() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.build();

		return confSpaces;
	}

	@Test
	public void test2RL0() {

		double epsilon = 0.95;
		Result result = runKStar(make2RL0(), epsilon, null, false, 1);
		assert2RL0(result, epsilon);
	}

	@Test
	public void test2RL0WithExternalMemory() {

		ExternalMemory.use(128, () -> {
			double epsilon = 0.95;
			Result result = runKStar(make2RL0(), epsilon, null, true, 1);
			assert2RL0(result, epsilon);
		});
	}

	private static void assert2RL0(Result result, double epsilon) {
		// check the results (values collected with e = 0.01 and 64 digits precision)
		// NOTE: these values don't match the ones in the TestKSImplLinear test because the conf spaces are slightly different
		// also, the new K* code has been updated to be more precise
		assertSequence(result,   0, "PHE ASP GLU THR PHE LYS ILE THR", 4.371921e+04, 4.470250e+30, 4.204899e+50, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [50.623756,50.627674] (log10)                    K* = 15.332751 in [15.332671,15.336670] (log10)
		assertSequence(result,   1, "PHE ASP GLU THR PHE LYS ILE SER", 4.371921e+04, 1.147556e+30, 4.049550e+50, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [30.059774,30.062460] (log10)                    complex [50.607407,50.611614] (log10)                    K* = 15.906961 in [15.904272,15.911168] (log10)
		assertSequence(result,   2, "PHE ASP GLU THR PHE LYS ILE ASN", 4.371921e+04, 4.842333e+29, 1.855842e+49, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [29.685055,29.685101] (log10)                    complex [49.268541,49.272112] (log10)                    K* = 14.942814 in [14.942765,14.946385] (log10)
		assertSequence(result,   3, "PHE ASP GLU THR PHE LYS ALA THR", 4.371921e+04, 1.600526e+27, 7.020108e+45, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [27.204263,27.204485] (log10)                    complex [45.846344,45.849857] (log10)                    K* = 14.001409 in [14.001185,14.004922] (log10)
		assertSequence(result,   4, "PHE ASP GLU THR PHE LYS VAL THR", 4.371921e+04, 5.970459e+28, 9.866975e+47, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [28.776008,28.776895] (log10)                    complex [47.994184,47.998118] (log10)                    K* = 14.577504 in [14.576615,14.581438] (log10)
		assertSequence(result,   5, "PHE ASP GLU THR PHE LYS LEU THR", 4.371921e+04, 3.784115e-11, 3.649388e+08, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [-10.42203,-10.42182] (log10)                    complex [8.562220 , 8.565874] (log10)                    K* = 14.343583 in [14.343374,14.347238] (log10)
		assertSequence(result,   6, "PHE ASP GLU THR PHE LYS PHE THR", 4.371921e+04, 2.881629e+24, null        , epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [24.459638,24.460021] (log10)                    complex [-Infinity,-Infinity] (log10,OutOfLowEnergies)   K* = none      in [-Infinity,-Infinity] (log10)
		assertSequence(result,   7, "PHE ASP GLU THR PHE LYS TYR THR", 4.371921e+04, 1.462260e+26, null        , epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [26.165025,26.167197] (log10)                    complex [-Infinity,-Infinity] (log10,OutOfLowEnergies)   K* = none      in [-Infinity,-Infinity] (log10)
		assertSequence(result,   8, "PHE ASP GLU THR PHE ASP ILE THR", 4.371921e+04, 4.358067e+20, 1.255853e+36, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [20.639294,20.639306] (log10)                    complex [36.098939,36.101958] (log10)                    K* = 10.818973 in [10.818958,10.821992] (log10)
		assertSequence(result,   9, "PHE ASP GLU THR PHE GLU ILE THR", 4.371921e+04, 4.898049e+20, 2.276285e+35, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [20.690023,20.690025] (log10)                    complex [35.357227,35.360904] (log10)                    K* = 10.026531 in [10.026527,10.030209] (log10)
		assertSequence(result,  10, "PHE ASP GLU THR TYR LYS ILE THR", 4.371921e+04, 4.703797e+30, 2.296368e+50, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [30.672449,30.673477] (log10)                    complex [50.361042,50.365003] (log10)                    K* = 15.047921 in [15.046890,15.051883] (log10)
		assertSequence(result,  11, "PHE ASP GLU THR ALA LYS ILE THR", 4.371921e+04, 3.343664e+28, 2.173203e+47, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [28.524223,28.524239] (log10)                    complex [47.337100,47.340468] (log10)                    K* = 14.172205 in [14.172187,14.175573] (log10)
		assertSequence(result,  12, "PHE ASP GLU THR VAL LYS ILE THR", 4.371921e+04, 9.149013e+29, 1.868156e+49, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [29.961374,29.961452] (log10)                    complex [49.271413,49.274723] (log10)                    K* = 14.669367 in [14.669286,14.672677] (log10)
		assertSequence(result,  13, "PHE ASP GLU THR ILE LYS ILE THR", 4.371921e+04, 3.475808e+30, 1.599260e+50, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [30.541056,30.541461] (log10)                    complex [50.203919,50.207685] (log10)                    K* = 15.022191 in [15.021784,15.025957] (log10)
		assertSequence(result,  14, "PHE ASP GLU THR LEU LYS ILE THR", 4.371921e+04, 5.495454e+27, 1.046077e+47, epsilon); // protein [4.640672 , 4.640674] (log10)                    ligand [27.740004,27.740644] (log10)                    complex [47.019564,47.023445] (log10)                    K* = 14.638888 in [14.638245,14.642769] (log10)
		assertSequence(result,  15, "PHE ASP GLU SER PHE LYS ILE THR", 3.283712e+06, 4.470250e+30, 1.478944e+53, epsilon); // protein [6.516365 , 6.518587] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [53.169952,53.174250] (log10)                    K* = 16.003255 in [16.000954,16.007553] (log10)
		assertSequence(result,  16, "PHE ASP GLU ASN PHE LYS ILE THR", 1.542601e+06, 4.470250e+30, 9.600623e+52, epsilon); // protein [6.188254 , 6.188610] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [52.982299,52.986487] (log10)                    K* = 16.143714 in [16.143280,16.147902] (log10)
		assertSequence(result,  17, "PHE ASP GLU GLN PHE LYS ILE THR", 2.630082e+06, 4.470250e+30, 3.349509e+53, epsilon); // protein [6.419969 , 6.420740] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [53.524981,53.529197] (log10)                    K* = 16.454680 in [16.453832,16.458896] (log10)
		assertSequence(result,  18, "PHE ASP ASP THR PHE LYS ILE THR", 1.227974e+01, 4.470250e+30, 1.473056e+45, epsilon); // protein [1.089189 , 1.089195] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [45.168219,45.171981] (log10)                    K* = 13.428698 in [13.428614,13.432460] (log10)
		assertSequence(result,  19, "PHE GLU GLU THR PHE LYS ILE THR", 2.026898e+05, 4.470250e+30, 1.098638e+50, epsilon); // protein [5.306832 , 5.306865] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [50.040855,50.044593] (log10)                    K* = 14.083691 in [14.083579,14.087429] (log10)
		assertSequence(result,  20, "TYR ASP GLU THR PHE LYS ILE THR", 1.699480e+04, 4.470250e+30, 2.827698e+46, epsilon); // protein [4.230316 , 4.230334] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [46.451433,46.455174] (log10)                    K* = 11.570785 in [11.570689,11.574526] (log10)
		assertSequence(result,  21, "ALA ASP GLU THR PHE LYS ILE THR", 6.128989e+02, 4.470250e+30, 1.679083e+45, epsilon); // protein [2.787389 , 2.787390] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [45.225072,45.228224] (log10)                    K* = 11.787351 in [11.787272,11.790503] (log10)
		assertSequence(result,  22, "VAL ASP GLU THR PHE LYS ILE THR", 1.273574e+02, 4.470250e+30, 2.389220e+45, epsilon); // protein [2.105024 , 2.105027] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [45.378256,45.381563] (log10)                    K* = 12.622900 in [12.622819,12.626207] (log10)
		assertSequence(result,  23, "ILE ASP GLU THR PHE LYS ILE THR", 6.030325e+02, 4.470250e+30, 2.015890e+46, epsilon); // protein [2.780341 , 2.780363] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [46.304467,46.308586] (log10)                    K* = 12.873794 in [12.873693,12.877914] (log10)
		assertSequence(result,  24, "LEU ASP GLU THR PHE LYS ILE THR", 4.638796e+00, 4.470250e+30, 4.750455e+43, epsilon); // protein [0.666405 , 0.666410] (log10)                    ligand [30.650332,30.650410] (log10)                    complex [43.676735,43.680566] (log10)                    K* = 12.359998 in [12.359915,12.363829] (log10)
	}

	public static ConfSpaces make1GUA11() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("1", "180")
			.build();
		protein.flexibility.get("21").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("24").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("25").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("27").setLibraryRotamers(Strand.WildType, "HID").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("29").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("40").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("181", "215")
			.build();
		ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("213").setLibraryRotamers(Strand.WildType, "HID", "HIE", "LYS", "ARG").addWildTypeRotamers().setContinuous();

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.build();

		return confSpaces;
	}

	@Test
	public void test1GUA11() {

		double epsilon = 0.999999;
		Result result = runKStar(make1GUA11(), epsilon, null, false, 1);

		// check the results (values collected with e = 0.1 and 64 digits precision)
		assertSequence(result,   0, "HIE VAL", 1.194026e+42, 2.932628e+07, 1.121625e+66, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [7.467257 , 7.467257] (log10)                    complex [66.049848,66.051195] (log10)                    K* = 16.505577 in [16.505576,16.506925] (log10)
		assertSequence(result,   1, "HIE HID", 1.194026e+42, 5.738568e+07, 3.346334e+66, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [7.758803 , 7.758803] (log10)                    complex [66.524569,66.543073] (log10)                    K* = 16.688752 in [16.688752,16.707256] (log10)
		assertSequence(result,   2, "HIE HIE", 1.194026e+42, 6.339230e+06, 5.544100e+65, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [6.802036 , 6.802036] (log10)                    complex [65.743831,65.769366] (log10)                    K* = 16.864781 in [16.864780,16.890316] (log10)
		assertSequence(result,   3, "HIE LYS", 1.194026e+42, 6.624443e+04, 3.315130e+63, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [4.821149 , 4.826752] (log10)                    complex [63.520501,63.563549] (log10)                    K* = 16.622337 in [16.616735,16.665386] (log10)
		assertSequence(result,   4, "HIE ARG", 1.194026e+42, 1.196619e+05, 5.375633e+64, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [5.077956 , 5.087238] (log10)                    complex [64.730430,64.774106] (log10)                    K* = 17.575460 in [17.566178,17.619136] (log10)
		assertSequence(result,   5, "HID VAL", 9.813429e+41, 2.932628e+07, 2.680104e+66, epsilon); // protein [41.991821,41.992159] (log10)                    ligand [7.467257 , 7.467257] (log10)                    complex [66.428152,66.446408] (log10)                    K* = 16.969074 in [16.968735,16.987330] (log10)
	}

	@Test
	public void test2RL0WithConfDB() {

		final double epsilon = 0.95;
		final String confdbPattern = "kstar.*.conf.db";
		final ConfSpaces confSpaces = make2RL0();

		try (TempFile proteinDBFile = new TempFile("kstar.protein.conf.db")) {
			try (TempFile ligandDBFile = new TempFile("kstar.ligand.conf.db")) {
				try (TempFile complexDBFile = new TempFile("kstar.complex.conf.db")) {

					// run with empty dbs
					Stopwatch sw = new Stopwatch().start();
					Result result = runKStar(confSpaces, epsilon, confdbPattern, false, 1);
					assert2RL0(result, epsilon);
					System.out.println(sw.getTime(2));

					// the dbs should have stuff in them

					try (ConfDB confdb = new ConfDB(confSpaces.protein, proteinDBFile)) {
						assertThat(confdb.getNumSequences(), greaterThan(0L));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					}

					try (ConfDB confdb = new ConfDB(confSpaces.ligand, ligandDBFile)) {
						assertThat(confdb.getNumSequences(), greaterThan(0L));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					}

					try (ConfDB confdb = new ConfDB(confSpaces.complex, complexDBFile)) {
						assertThat(confdb.getNumSequences(), greaterThan(0L));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					}

					assertThat(proteinDBFile.exists(), is(true));
					assertThat(ligandDBFile.exists(), is(true));
					assertThat(complexDBFile.exists(), is(true));

					// run again with full dbs
					sw = new Stopwatch().start();
					Result result2 = runKStar(confSpaces, epsilon, confdbPattern, false, 1);
					assert2RL0(result2, epsilon);
					System.out.println(sw.getTime(2));
				}
			}
		}
	}

	public static ConfSpaces make2RL0OnlyOneMutant() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

		// define the protein strand with just wild-type
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType).setContinuous();

		// define the ligand strand with one mutant, and no wild-type
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A193").setLibraryRotamers("VAL").setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.build();

		return confSpaces;
	}

	@Test
	public void test2RL0OnlyOneMutant() {

		double epsilon = 0.99;
		Result result = runKStar(make2RL0OnlyOneMutant(), epsilon, null, false, 1);

		assertThat(result.scores.size(), is(1));
		assertThat(result.scores.get(0).sequence.toString(Sequence.Renderer.AssignmentMutations), is("A193=VAL"));
	}

	public static ConfSpaces make2RL0SpaceWithoutWildType() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

		// define the protein strand with just wild-type
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "VAL").setContinuous();

		// define the ligand strand with one mutant, and no wild-type
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A193").setLibraryRotamers("VAL").setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.build();

		return confSpaces;
	}

	@Test
	public void test2RL0SpaceWithoutWildType() {

		double epsilon = 0.99;
		Result result = runKStar(make2RL0SpaceWithoutWildType(), epsilon, null, false, 2);

		assertThat(result.scores.size(), is(2));
		assertThat(result.scores.get(0).sequence.toString(Sequence.Renderer.AssignmentMutations), is("G654=VAL A193=VAL"));
		assertThat(result.scores.get(1).sequence.toString(Sequence.Renderer.AssignmentMutations), is("G654=thr A193=VAL"));
	}

	public static void assertSequence(Result result, int sequenceIndex, String sequence, Double proteinQStar, Double ligandQStar, Double complexQStar, double epsilon) {

		KStar.ScoredSequence scoredSequence = result.scores.get(sequenceIndex);

		// check the sequence
		assertThat(scoredSequence.sequence.toString(Sequence.Renderer.ResType), is(sequence));

		// check q* values and epsilon
		assertResult(scoredSequence.score.protein, proteinQStar, epsilon);
		assertResult(scoredSequence.score.ligand, ligandQStar, epsilon);
		assertResult(scoredSequence.score.complex, complexQStar, epsilon);
	}

	public static void assertResult(PartitionFunction.Result result, Double qstar, double epsilon) {
		if (qstar != null) {
			assertThat(result.status, is(PartitionFunction.Status.Estimated));
			assertThat(result.values.qstar.doubleValue(), greaterThanOrEqualTo(qstar*(1.0 - epsilon)));
			assertThat(result.values.getEffectiveEpsilon(), lessThanOrEqualTo(epsilon));
		} else {
			assertThat(result.status, is(not(PartitionFunction.Status.Estimated)));
		}
	}
}
