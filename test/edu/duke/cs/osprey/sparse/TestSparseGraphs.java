package edu.duke.cs.osprey.sparse;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;

public class TestSparseGraphs {


	private Problem problemDiscrete;
	private Problem problemContinuous;
	private Problem problemMultipleStrands;
	private Problem problemBigContinuous;


	private static class Problem {

		public final SimpleConfSpace confSpace;
		public final ForcefieldParams ffparams;
		public final EnergyCalculator ecalc;
		public final ConfEnergyCalculator confEcalc;
		public final EnergyMatrix emat;

		public Problem(Strand ... strands) {
			confSpace = new SimpleConfSpace.Builder().addStrands(strands).build();
			ffparams = new ForcefieldParams();
			ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build();
			confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}

		public SimpleGMECFinder makeFinder() {
			return new SimpleGMECFinder.Builder(
				new ConfAStarTree.Builder(emat, confSpace).build(),
				confEcalc
			)
			.build();
		}

		public SimpleGMECFinder makeExternalFinder() {
			return new SimpleGMECFinder.Builder(
				new ConfAStarTree.Builder(emat, confSpace)
					.useExternalMemory()
					.build(),
				confEcalc
			)
			.useExternalMemory()
			.build();
		}

		public SimpleGMECFinder makeConfDBFinder(File confdbFile, Integer interruptAtConfNum) {
			return new SimpleGMECFinder.Builder(
				new ConfAStarTree.Builder(emat, confSpace).build(),
				new ConfEnergyCalculator(confEcalc) {

					int numConfs = 0;

					// override the energy calculation to throw an error after a few confs
					@Override
					public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag, ResidueInteractions inters) {

						if (interruptAtConfNum != null && numConfs++ == interruptAtConfNum) {
							throw new Error("Interrupted!");
						}

						return confEcalc.calcEnergy(frag, inters);
					}
				}
			)
			.setConfDB(confdbFile)
			.build();
		}
	}
    @Test
    public void testGenerateGraphs() {
        System.out.println("Test!");

		Molecule mol = PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb");

		Strand strand1 = new Strand.Builder(mol).build();
		strand1.flexibility.get("A2").setLibraryRotamers("ALA", "GLY");
		strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "VAL");
		strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		problemDiscrete = new Problem(strand1);

		Strand strand2 = new Strand.Builder(mol).setResidues("A2", "A30").build();
		strand2.flexibility.get("A2").setLibraryRotamers("ALA", "GLY");
		strand2.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "VAL", "ARG").setContinuous(10);
		strand2.flexibility.get("A4").addWildTypeRotamers();
		problemContinuous = new Problem(strand2);

		Strand strand3 = new Strand.Builder(mol).setResidues("A2", "A30").build();
		strand3.flexibility.get("A2").addWildTypeRotamers();
		strand3.flexibility.get("A3").addWildTypeRotamers();
		strand3.flexibility.get("A4").addWildTypeRotamers();
		Strand strand4 = new Strand.Builder(mol).setResidues("A31", "A60").build();
		strand4.flexibility.get("A31").addWildTypeRotamers();
		strand4.flexibility.get("A32").addWildTypeRotamers();
		strand4.flexibility.get("A33").addWildTypeRotamers();
		problemMultipleStrands = new Problem(strand3, strand4);

		Strand strand5 = new Strand.Builder(mol).setResidues("A2", "A30").build();
		for (String resNum : Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7")) {
			strand5.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType, "ALA", "GLY")
				.setContinuous();
		}
		problemBigContinuous = new Problem(strand5);
	}
}
