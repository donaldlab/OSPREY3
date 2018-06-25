package edu.duke.cs.osprey.newEwakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.ewakstar.EwakstarDoer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.Set;

import static edu.duke.cs.osprey.tools.Log.log;


public class EwakstarLab {

	public static void main(String[] args) {
		// run COMETS
		//EwakstarDoer ewakstarDoer = run2RL0();
		EwakstarDoer ewakstarDoer = run1GWC();
		Set<Sequence> seqs = ewakstarDoer.run(ewakstarDoer.state);
	}

	public static EwakstarDoer run2RL0(){

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		// especially since all the state conf spaces will add/share wild-type rotamers to/from the library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("G648", "G654")
				.build();
		ligand.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "ILE","VAL","LEU","ALA","TYR").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("G654").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A155", "A194")
				.build();
		protein.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the COMETS states

		EwakstarDoer.State PL = new EwakstarDoer.State(
				"PL",
				new SimpleConfSpace.Builder()
						.addStrand(protein)
						.addStrand(ligand)
						.build()
		);

		int orderMag = 10; //order of magnitude worse in partition function we want to keep sequences relative to the wild-type
		int numEWAKStarSeqs = 10000; //number of sequences we want to limit ourselves to during the "sequence filter" portion of ewakstarDoer
		double ewakstarEw = 30.0; //energy window within the wild-type for finding sequences in the "sequence filter" portion of ewakstarDoer
		int numPfConfs = 5000; //num of conformations for the partition function calculation
		double pfEw = 1.0; //partition function energy window calculation
		int numTopOverallSeqs = 5; //end result number of sequences we want K* estimates for
		int numCpus = 4;
		double epsilon = 0.01;

		EwakstarDoer ewakstarDoer = new EwakstarDoer.Builder()
				.setNumEWAKStarSeqs(numEWAKStarSeqs)
				.setOrderOfMag(orderMag)
				.setPfEw(pfEw)
				.setEpsilon(epsilon)
				.setNumPfConfs(numPfConfs)
				.setNumTopOverallSeqs(numTopOverallSeqs)
				.addState(PL)
				.setEw(ewakstarEw)
				.setMutableType("exact")
				.setNumMutable(1)
				.setSeqFilterOnly(false)
				.setNumCpus(numCpus)
				.setLogFile(new File("ewakstar.sequences.tsv"))
				.build();


		EnergyCalculator ecalc = new EnergyCalculator.Builder(PL.confSpace, ffparams)
				.setParallelism(Parallelism.makeCpu(numCpus))
				.build();
		EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(PL.confSpace, ffparams)
				.setIsMinimizing(false)
				.setParallelism(Parallelism.makeCpu(4))
				.build();

		// what are conformation energies?
		SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(PL.confSpace, ecalc)
				.build()
				.calcReferenceEnergies();
		SimpleReferenceEnergies rigidEref = new SimplerEnergyMatrixCalculator.Builder(PL.confSpace, rigidEcalc)
				.build()
				.calcReferenceEnergies();

		PL.confEcalc =new ConfEnergyCalculator.Builder(PL.confSpace,ecalc)
				.setReferenceEnergies(eref)
				.build();

		PL.confRigidEcalc =new ConfEnergyCalculator.Builder(PL.confSpace,rigidEcalc)
				.setReferenceEnergies(rigidEref)
				.build();

		// calc the energy matrix
		PL.emat =new SimplerEnergyMatrixCalculator.Builder(PL.confEcalc)
				.setCacheFile(new File(String.format("ewakstar.%s.emat", PL.name)))
				.build()
				.calcEnergyMatrix();

		PL.fragmentEnergies =PL.emat;
		PL.ematRigid =new SimplerEnergyMatrixCalculator.Builder(PL.confRigidEcalc)
				.setCacheFile(new File(String.format("ewakstar.%s.ematRigid", PL.name)))
				.build()
				.calcEnergyMatrix();

		// make the conf tree factory
		PL.confTreeFactory =(rcs)->new ConfAStarTree.Builder(PL.emat,rcs)
				.setTraditional()
				.build();

		return ewakstarDoer;
	}

	public static EwakstarDoer run1GWC(){

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/1GWC.shell.pdb");

		// make sure all strands share the same template library
		// especially since all the state conf spaces will add/share wild-type rotamers to/from the library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("B89", "B161")
				.build();
		ligand.flexibility.get("B93").setLibraryRotamers("TYR").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B100").setLibraryRotamers("PHE").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B101").setLibraryRotamers("TRP").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B104").setLibraryRotamers("TYR").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("C7", "C80")
				.build();
		protein.flexibility.get("C53").setLibraryRotamers("HIE", "ILE").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("C64").setLibraryRotamers("ALA").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("C67").setLibraryRotamers("SER", "CYS").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("C71").setLibraryRotamers("ILE").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("C78").setLibraryRotamers("GLU").addWildTypeRotamers().setContinuous();

		// make the COMETS states

		EwakstarDoer.State PL = new EwakstarDoer.State(
				"PL",
				new SimpleConfSpace.Builder()
						.addStrand(protein)
						.addStrand(ligand)
						.build()
		);

		int orderMag = 5; //order of magnitude worse in partition function we want to keep sequences relative to the wild-type
		int numEWAKStarSeqs = 10000; //number of sequences we want to limit ourselves to during the "sequence filter" portion of ewakstarDoer
		double ewakstarEw = 10.0; //energy window within the wild-type for finding sequences in the "sequence filter" portion of ewakstarDoer
		int numPfConfs = 500; //num of conformations for the partition function calculation
		double pfEw = 1.0; //partition function energy window calculation
		int numTopOverallSeqs = 5; //end result number of sequences we want K* estimates for
		int numCpus = 4;
		double epsilon = 0.68;

		EwakstarDoer ewakstarDoer = new EwakstarDoer.Builder()
				.setNumEWAKStarSeqs(numEWAKStarSeqs)
				.setOrderOfMag(orderMag)
				.setPfEw(pfEw)
				.setEpsilon(epsilon)
				.setNumPfConfs(numPfConfs)
				.setNumTopOverallSeqs(numTopOverallSeqs)
				.addState(PL)
				.setEw(ewakstarEw)
				.setMutableType("exact")
				.setNumMutable(1)
				.setSeqFilterOnly(false)
				.setNumCpus(numCpus)
				.setLogFile(new File("ewakstar.sequences.tsv"))
				.build();


		EnergyCalculator ecalc = new EnergyCalculator.Builder(PL.confSpace, ffparams)
				.setParallelism(Parallelism.makeCpu(numCpus))
				.build();
		EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(PL.confSpace, ffparams)
				.setIsMinimizing(false)
				.setParallelism(Parallelism.makeCpu(4))
				.build();

		// what are conformation energies?
		SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(PL.confSpace, ecalc)
				.build()
				.calcReferenceEnergies();
		SimpleReferenceEnergies rigidEref = new SimplerEnergyMatrixCalculator.Builder(PL.confSpace, rigidEcalc)
				.build()
				.calcReferenceEnergies();

		PL.confEcalc =new ConfEnergyCalculator.Builder(PL.confSpace,ecalc)
				.setReferenceEnergies(eref)
				.build();

		PL.confRigidEcalc =new ConfEnergyCalculator.Builder(PL.confSpace,rigidEcalc)
				.setReferenceEnergies(rigidEref)
				.build();

		// calc the energy matrix
		PL.emat =new SimplerEnergyMatrixCalculator.Builder(PL.confEcalc)
				.setCacheFile(new File(String.format("ewakstar.%s.emat", PL.name)))
				.build()
				.calcEnergyMatrix();

		PL.fragmentEnergies =PL.emat;
		PL.ematRigid =new SimplerEnergyMatrixCalculator.Builder(PL.confRigidEcalc)
				.setCacheFile(new File(String.format("ewakstar.%s.ematRigid", PL.name)))
				.build()
				.calcEnergyMatrix();


		// make the conf tree factory
		PL.confTreeFactory =(rcs)->new ConfAStarTree.Builder(PL.emat,rcs)
				.setTraditional()
				.build();

		return ewakstarDoer;
	}

}
