package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.MinimizingConfEnergyCalculator;
import edu.duke.cs.osprey.energy.MinimizingFragmentEnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;

public class EnergyPartitionsPlayground {
	
	public static void main(String[] args) {
		
		// read a protein
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb")).build();
		
		// configure flexibility
		/*
		for (int i=2; i<=15; i++) {
			strand.flexibility.get(i).setLibraryRotamers("LEU", "ILE").setContinuous();
		}
		*/
		for (int i=2; i<=8; i++) {
			strand.flexibility.get(i).setLibraryRotamers("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "CYS", "ASN", "GLN", "ASP", "GLU", "PHE", "TRP", "TYR", "HIE", "HID", "LYS", "ARG", "MET").setContinuous();
		}
		
		// make the conf space
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand)
			//.setShellDistance(9)
			.build();
		
		System.out.println("Residues: " + strand.mol.residues.size());
		System.out.println("Shell residues: " + confSpace.shellResNumbers.size());
		System.out.println("positions: " + confSpace.positions.size());
		System.out.println("RCs per pos: " + confSpace.positions.get(0).resConfs.size());
		
		// choose the default forcefield
		ForcefieldParams ffparams = new ForcefieldParams();
		
		MinimizingFragmentEnergyCalculator fragEcalc = new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeGpu(1, 16))
			.setType(MinimizingFragmentEnergyCalculator.Type.ResidueCudaCCD)
			.build();
		
		// get reference energies
		SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, fragEcalc)
			.build()
			.calcReferenceEnergies();

		// compute the energy matrix
		EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, fragEcalc)
			.setReferenceEnergies(eref)
			//.setEnergyPartition(EnergyPartition.Traditional)
			.setEnergyPartition(EnergyPartition.AllOnPairs)
			.build()
			.calcEnergyMatrix();
		
		// config TPIE
		ExternalMemory.setInternalLimit(128);
		ExternalMemory.setTempDir(System.getProperty("user.home"), "tpie-epart");
		
		// how should confs be ordered?
		ConfSearch confSearch = new ConfAStarTree.Builder(emat, confSpace)
			//.setTraditional()
			.setMPLP(new ConfAStarTree.MPLPBuilder().setUpdater(new EdgeUpdater()).setNumIterations(5))
			.setUseExternalMemory()
			.setShowProgress(true)
			.build();
	
		// what's the energy of a conformation?
		MinimizingConfEnergyCalculator confEcalc = new MinimizingConfEnergyCalculator.Builder(fragEcalc)
			.setReferenceEnergies(eref)
			.build();
		
		// find the GMEC!
		System.out.println("Finding GMEC...");
		EnergiedConf gmec = new SimpleGMECFinder.Builder(confSpace, confSearch, confEcalc)
			.setPrintIntermediateConfsToConsole(false)
			.setUseExternalMemory()
			.build()
			.find();
		
		/* for the original design
		// check the GMEC is still correct
		//assertThat(gmec.getAssignments(), is(new int[] { 6, 3, 9, 4, 1 }));
		//assertThat(gmec.getEnergy(), isAbsolutely(-20.891946, 1e-6));
		
		// GMEC score:      -21.623094 (gap: ~0.73)
		// min-bound score: -22.305976 (gap: ~8.89) !!!
		*/
	}
}
