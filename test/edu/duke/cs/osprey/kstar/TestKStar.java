package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import java.util.concurrent.atomic.AtomicReference;

public class TestKStar {

	// TODO: get rid of eref setting
	public static KStar runKStar(ForcefieldParams ffparams, SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, double epsilon) {

		AtomicReference<KStar> kstarRef = new AtomicReference<>(null);

		// TEMP
		//Parallelism parallelism = Parallelism.makeCpu(4);
		Parallelism parallelism = Parallelism.make(4, 1, 8);

		// how should we compute energies of molecules?
		new EnergyCalculator.Builder(complex, ffparams)
			.setParallelism(parallelism)
			.use((ecalc) -> {

				// how should we define energies of conformations?
				KStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
					return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
						.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
							.build()
							.calcReferenceEnergies()
						).build();
				};

				// how should confs be ordered and searched?
				ConfSearchFactory confSearchFactory = (emat, pmat) -> {
					return new ConfAStarTree.Builder(emat, pmat)
						.setTraditional()
						.build();
				};

				// run K*
				KStar.Settings settings = new KStar.Settings.Builder()
					.setEpsilon(epsilon)
					.addScoreConsoleWriter(new KStarScoreWriter.Formatter.Test())
					//.setShowPfuncProgress(true)
					.build();
				KStar kstar = new KStar(protein, ligand, complex, ecalc, confEcalcFactory, confSearchFactory, settings);
				kstar.run();

				// pass back the ref
				kstarRef.set(kstar);
			});

		return kstarRef.get();
	}

	@Test
	public void test2RL0() {

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(648, 654)
			.build();
		protein.flexibility.get(649).setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get(650).setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get(651).setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		protein.flexibility.get(654).setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(155, 194)
			.build();
		ligand.flexibility.get(156).setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(172).setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(192).setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(193).setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		SimpleConfSpace proteinConfSpace = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		SimpleConfSpace ligandConfSpace = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();
		SimpleConfSpace complexConfSpace = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.build();

		double epsilon = 0.95;
		KStar kstar = runKStar(ffparams, proteinConfSpace, ligandConfSpace, complexConfSpace, epsilon);

		// check the results (values collected with e = 0.01)
		// NOTE: these values don't match the ones in the TestKSImplLinear test because the conf spaces are slightly different
		assertSequence(kstar,   0, "PHE ASP GLU THR PHE LYS ILE THR", 2.797907e+04, 4.445880e+30, 4.203126e+50, epsilon);
		assertSequence(kstar,   1, "PHE ASP GLU THR PHE LYS ILE SER", 2.797907e+04, 1.145512e+30, 4.047880e+50, epsilon);
		assertSequence(kstar,   2, "PHE ASP GLU THR PHE LYS ILE ASN", 2.797907e+04, 4.826033e+29, 1.855048e+49, epsilon);
		assertSequence(kstar,   3, "PHE ASP GLU THR PHE LYS ALA THR", 2.797907e+04, 1.592319e+27, 7.016956e+45, epsilon);
		assertSequence(kstar,   4, "PHE ASP GLU THR PHE LYS VAL THR", 2.797907e+04, 5.941138e+28, 9.862751e+47, epsilon);
		assertSequence(kstar,   5, "PHE ASP GLU THR PHE LYS LEU THR", 2.797907e+04, null,         3.647796e+08, epsilon);
		assertSequence(kstar,   6, "PHE ASP GLU THR PHE LYS PHE THR", 2.797907e+04, 2.871398e+24, null,         epsilon);
		assertSequence(kstar,   7, "PHE ASP GLU THR PHE LYS TYR THR", 2.797907e+04, 1.459454e+26, null,         epsilon);
		assertSequence(kstar,   8, "PHE ASP GLU THR PHE ASP ILE THR", 2.797907e+04, 4.343893e+20, 1.255309e+36, epsilon);
		assertSequence(kstar,   9, "PHE ASP GLU THR PHE GLU ILE THR", 2.797907e+04, 4.885325e+20, 2.275234e+35, epsilon);
		assertSequence(kstar,  10, "PHE ASP GLU THR TYR LYS ILE THR", 2.797907e+04, 4.686693e+30, 2.295400e+50, epsilon);
		assertSequence(kstar,  11, "PHE ASP GLU THR ALA LYS ILE THR", 2.797907e+04, 3.331113e+28, 2.172265e+47, epsilon);
		assertSequence(kstar,  12, "PHE ASP GLU THR VAL LYS ILE THR", 2.797907e+04, 9.112545e+29, 1.867349e+49, epsilon);
		assertSequence(kstar,  13, "PHE ASP GLU THR ILE LYS ILE THR", 2.797907e+04, 3.460857e+30, 1.598583e+50, epsilon);
		assertSequence(kstar,  14, "PHE ASP GLU THR LEU LYS ILE THR", 2.797907e+04, 5.472100e+27, 1.045633e+47, epsilon);
		assertSequence(kstar,  15, "PHE ASP GLU SER PHE LYS ILE THR", 2.108199e+06, 4.445880e+30, 1.478335e+53, epsilon);
		assertSequence(kstar,  16, "PHE ASP GLU ASN PHE LYS ILE THR", 9.882385e+05, 4.445880e+30, 9.596655e+52, epsilon);
		assertSequence(kstar,  17, "PHE ASP GLU GLN PHE LYS ILE THR", 1.686659e+06, 4.445880e+30, 3.348091e+53, epsilon);
		assertSequence(kstar,  18, "PHE ASP ASP THR PHE LYS ILE THR", 6.687589e+00, 4.445880e+30, 1.472011e+45, epsilon);
		assertSequence(kstar,  19, "PHE GLU GLU THR PHE LYS ILE THR", 1.283597e+05, 4.445880e+30, 1.098175e+50, epsilon);
		assertSequence(kstar,  20, "TYR ASP GLU THR PHE LYS ILE THR", 1.698008e+04, 4.445880e+30, 2.827534e+46, epsilon);
		assertSequence(kstar,  21, "ALA ASP GLU THR PHE LYS ILE THR", 6.128158e+02, 4.445880e+30, 1.678916e+45, epsilon);
		assertSequence(kstar,  22, "VAL ASP GLU THR PHE LYS ILE THR", 1.273320e+02, 4.445880e+30, 2.389072e+45, epsilon);
		assertSequence(kstar,  23, "ILE ASP GLU THR PHE LYS ILE THR", 6.026501e+02, 4.445880e+30, 2.015879e+46, epsilon);
		assertSequence(kstar,  24, "LEU ASP GLU THR PHE LYS ILE THR", 4.636954e+00, 4.445880e+30, 4.750328e+43, epsilon);
	}

	@Test
	public void test1GUA_11() {

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(1, 180)
			.build();
		protein.flexibility.get(21).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get(24).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get(25).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get(27).setLibraryRotamers(Strand.WildType, "HID").addWildTypeRotamers().setContinuous();
		protein.flexibility.get(29).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get(40).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(181, 215)
			.build();
		ligand.flexibility.get(209).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(213).setLibraryRotamers(Strand.WildType, "HID", "HIE", "LYS", "ARG").addWildTypeRotamers().setContinuous();

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		double shellDist = 10.0;
		SimpleConfSpace proteinConfSpace = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.setShellDistance(shellDist)
			.build();
		SimpleConfSpace ligandConfSpace = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.setShellDistance(shellDist)
			.build();
		SimpleConfSpace complexConfSpace = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.setShellDistance(shellDist)
			.build();

		// TEMP
		//double epsilon = 0.99;
		double epsilon = 0.4;
		KStar kstar = runKStar(ffparams, proteinConfSpace, ligandConfSpace, complexConfSpace, epsilon);

		// TODO NEXTTIME: doesn't quite match these results yet
		// conf space doesn't seem to be the problem... must be a setting somewhere?
		// need to check pfunc conf scores and energies, start with wildtype I guess...

		// check the results (values collected with e = 0.1)
		assertSequence(kstar,   0, "ILE ILE GLN HIE VAL TYR LYS VAL", 1.1838e+42, 2.7098e+7, 1.1195e+66, epsilon);
		assertSequence(kstar,   1, "ILE ILE GLN HIE VAL TYR LYS HID", 1.1838e+42, 5.5334e+7, 3.3455e+66, epsilon);
		assertSequence(kstar,   2, "ILE ILE GLN HIE VAL TYR LYS HIE", 1.1838e+42, 5.8485e+6, 5.5426e+65, epsilon);
		assertSequence(kstar,   3, "ILE ILE GLN HIE VAL TYR LYS LYS", 1.1838e+42, 6.3856e+4, 3.3162e+63, epsilon);
		assertSequence(kstar,   4, "ILE ILE GLN HIE VAL TYR LYS ARG", 1.1838e+42, 1.1527e+5, 5.3772e+64, epsilon);
		assertSequence(kstar,   5, "ILE ILE GLN HID VAL TYR LYS VAL", 9.7448e+41, 2.7098e+7, 2.6775e+66, epsilon);
	}

	public static void assertSequence(KStar kstar, int sequenceIndex, String sequence, Double proteinQStar, Double ligandQStar, Double complexQStar, double epsilon) {

		// check the sequence
		assertThat(kstar.complex.sequences.get(sequenceIndex), is(new KStar.Sequence(sequence)));

		// check q* values and epsilon
		assertResult(kstar.protein.pfuncResults.get(kstar.protein.sequences.get(sequenceIndex)), proteinQStar, epsilon);
		assertResult(kstar.ligand.pfuncResults.get(kstar.ligand.sequences.get(sequenceIndex)), ligandQStar, epsilon);
		assertResult(kstar.complex.pfuncResults.get(kstar.complex.sequences.get(sequenceIndex)), complexQStar, epsilon);
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
