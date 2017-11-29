package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import java.io.File;
import java.util.Random;
import static edu.duke.cs.osprey.TestBase.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.kstar.TestKStar.ConfSpaces;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

public class TestSequenceEnsembleAnalyzer {

	private static final double EnergyEpsilon = 1e-6;
	
	private static TestKStar.ConfSpaces makeConfSpaces() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/gp120SRVRC26.09SR.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(mol)
			.addTemplates(FileTools.readFile("examples/python.KStar/all_nuc94_and_gr.tys.in"))
			.addTemplateCoords(FileTools.readFile("examples/python.KStar/all_amino_coords.tys.in"))
			.addRotamers(FileTools.readFile("examples/python.KStar/GenericRotamers.tys.dat"))
			.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("H1792", "L2250")
			.build();

		protein.flexibility.get("H1901").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("H1904").setLibraryRotamers("ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", 
				"SER", "THR", "LYS", "ARG", "HIS", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("H1905").setLibraryRotamers("ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET",
				"SER", "THR", "LYS", "ARG", "HIS", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("H1906").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("H1907").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("H1908").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("F379", "J1791")
			.build();

		ligand.flexibility.get("G973").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("G977").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("G978").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("G979").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("G980").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("J1448").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
			.addStrands(ligand, protein)
			.build();

		return confSpaces;
	}

	@Test
	public void testHIV() {

		TestKStar.ConfSpaces confSpaces = makeConfSpaces();

		Parallelism parallelism = Parallelism.makeCpu(4);
		//Parallelism parallelism = Parallelism.make(4, 1, 8);
		ExternalMemory.setInternalLimit(500);

		// how should we compute energies of molecules?
		new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
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
				KStar.ConfSearchFactory confSearchFactory = (emat, pmat) -> {
					return new ConfAStarTree.Builder(emat, pmat)
						.useExternalMemory()
						.setTraditional()
						.build();
				};

				KStar.Settings settings = new KStar.Settings.Builder()
						.setEpsilon(0.68)
						.setEnergyMatrixCachePattern("*.emat")
						.build();

				SequenceEnsembleAnalyzer analyzer = new SequenceEnsembleAnalyzer(
					confSpaces.protein,
					confSpaces.ligand,
					confSpaces.complex,
					ecalc,
					confEcalcFactory,
					confSearchFactory,
					settings
				);

				analyzer.analyze(
					new KStar.Sequence("THR ARG ASP LYS LYS ARG ASP ASP TYR GLY LYS GLN"),
					10000
				);


			});
	}
}
