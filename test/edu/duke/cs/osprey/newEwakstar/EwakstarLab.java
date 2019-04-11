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
import edu.duke.cs.osprey.pruning.PruningMatrix;
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
		EwakstarDoer ewakstarDoer = run2RL0();
		//EwakstarDoer ewakstarDoer = run1GUA();
		//EwakstarDoer ewakstarDoer = run1GWC();
		//EwakstarDoer ewakstarDoer = runSpA();
		//EwakstarDoer ewakstarDoer = run2HNV();
		Set<Sequence> seqs = ewakstarDoer.run(ewakstarDoer.state);
	}

	public static EwakstarDoer run2HNV(){

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2hnv_prepped.pdb");

		// make sure all strands share the same template library
		// especially since all the state conf spaces will add/share wild-type rotamers to/from the library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A7", "A87")
				.build();
		ligand.flexibility.get("A25").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIP", "HIE", "HID", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A36").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("B7", "B87")
				.build();
		protein.flexibility.get("B81").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("B87").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the COMETS states

		EwakstarDoer.State PL = new EwakstarDoer.State(
				"PL",
				new SimpleConfSpace.Builder()
						.addStrand(protein)
						.addStrand(ligand)
						.build()
		);

		int orderMag = 3; //order of magnitude worse in partition function we want to keep sequences relative to the wild-type
		double ewakstarEw = 30.0; //energy window within the wild-type for finding sequences in the "sequence filter" portion of ewakstarDoer
		int numPfConfs = 5000; //num of conformations for the partition function calculation
		double pfEw = 1.0; //partition function energy window calculation
		int numTopOverallSeqs = 5; //end result number of sequences we want K* estimates for
		int numCpus = 4;
		double epsilon = 0.1;

		EwakstarDoer ewakstarDoer = new EwakstarDoer.Builder()
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

		PruningMatrix pmat = new SimpleDEE.Runner().run(PL.confSpace,PL.emat);

		PL.fragmentEnergies =PL.emat;
		PL.ematRigid =new SimplerEnergyMatrixCalculator.Builder(PL.confRigidEcalc)
				.setCacheFile(new File(String.format("ewakstar.%s.ematRigid", PL.name)))
				.build()
				.calcEnergyMatrix();

		// make the conf tree factory
		PL.confTreeFactoryMin =(rcs)->new ConfAStarTree.Builder(PL.emat,rcs)
				.setMaxNumNodes(20000000)
				.setTraditional()
				.build();

		PL.confTreeFactoryRigid =(rcs)->new ConfAStarTree.Builder(PL.ematRigid,rcs)
				.setMaxNumNodes(20000000)
				.setTraditional()
				.build();

		return ewakstarDoer;
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

		int orderMag = 3; //order of magnitude worse in partition function we want to keep sequences relative to the wild-type
		double ewakstarEw = 30.0; //energy window within the wild-type for finding sequences in the "sequence filter" portion of ewakstarDoer
		int numPfConfs = 5000; //num of conformations for the partition function calculation
		double pfEw = 1.0; //partition function energy window calculation
		int numTopOverallSeqs = 5; //end result number of sequences we want K* estimates for
		int numCpus = 4;
		double epsilon = 0.01;

		EwakstarDoer ewakstarDoer = new EwakstarDoer.Builder()
				.setOrderOfMag(orderMag)
				.setPrintPDBs(false)
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
		PL.confTreeFactoryMin =(rcs)->new ConfAStarTree.Builder(PL.emat,rcs)
				.setMaxNumNodes(20000000)
				.setTraditional()
				.build();

		PL.confTreeFactoryRigid =(rcs)->new ConfAStarTree.Builder(PL.ematRigid,rcs)
				.setMaxNumNodes(20000000)
				.setTraditional()
				.build();

		return ewakstarDoer;
	}

	public static EwakstarDoer run1GUA(){

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/1gua_adj.min.pdb");

		// make sure all strands share the same template library
		// especially since all the state conf spaces will add/share wild-type rotamers to/from the library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("181", "215")
				.build();
		ligand.flexibility.get("182").setLibraryRotamers(Strand.WildType, "LYS").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("184").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("194").setLibraryRotamers(Strand.WildType, "LYS").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("196").setLibraryRotamers(Strand.WildType, "ASN", "VAL").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("208").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("210").setLibraryRotamers(Strand.WildType, "SER", "LYS").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("212").setLibraryRotamers(Strand.WildType, "TYR").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("213").setLibraryRotamers(Strand.WildType, "TYR").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("1", "180")
				.build();
		protein.flexibility.get("21").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("24").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("27").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("29").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("36").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("37").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("38").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("40").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the COMETS states

		EwakstarDoer.State PL = new EwakstarDoer.State(
				"PL",
				new SimpleConfSpace.Builder()
						.addStrand(protein)
						.addStrand(ligand)
						.build()
		);

		boolean useSMA = false;
		int orderMag = 3; //order of magnitude worse in partition function we want to keep sequences relative to the wild-type
		double ewakstarEw = 2.0; //energy window within the wild-type for finding sequences in the "sequence filter" portion of ewakstarDoer
		int numPfConfs = 1; //num of conformations for the partition function calculation
		double pfEw = 1.0; //partition function energy window calculation
		int numTopOverallSeqs = 5; //end result number of sequences we want K* estimates for
		int numCpus = 1;
		double epsilon = 0.99;

		EwakstarDoer ewakstarDoer = new EwakstarDoer.Builder()
				.setPrintPDBs(true)
				.setUseWtBenchmark(false)
				.setOrderOfMag(orderMag)
				.setPfEw(pfEw)
				.setEpsilon(epsilon)
				.setNumPfConfs(numPfConfs)
				.setNumTopOverallSeqs(numTopOverallSeqs)
				.addState(PL)
				.setEw(ewakstarEw)
				.setMutableType("exact")
				.setNumMutable(1)
				.setSeqFilterOnly(true)
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

		if (useSMA) {
			PL.confTreeFactoryMin = (rcs) -> new ConfAStarTree.Builder(PL.emat, rcs)
					.setMaxNumNodes(20000000)
					.setTraditional()
					.build();

			PL.confTreeFactoryRigid = (rcs) -> new ConfAStarTree.Builder(PL.ematRigid, rcs)
					.setMaxNumNodes(20000000)
					.setTraditional()
					.build();
		} else{
			PL.confTreeFactoryMin = (rcs) -> new ConfAStarTree.Builder(PL.emat, rcs)
					.setTraditional()
					.build();

			PL.confTreeFactoryRigid = (rcs) -> new ConfAStarTree.Builder(PL.ematRigid, rcs)
					.setTraditional()
					.build();
		}

		return ewakstarDoer;
	}

	public static EwakstarDoer runSpA(){
		Molecule mol = PDBIO.readFile("./test-resources/4npd.A_DomainA_noH_trim_his_clean.min.pdb");
		ForcefieldParams ffparams = new ForcefieldParams();
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		Strand protein = new Strand.Builder(mol)
				.setResidues("A1", "A36")
				.setTemplateLibrary(templateLib)
				.build();
		protein.flexibility.get("A8").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


		Strand ligand = new Strand.Builder(mol)
				.setResidues("A37", "A58")
				.setTemplateLibrary(templateLib)
				.build();
		ligand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A42").setLibraryRotamers(Strand.WildType, "THR", "LYS").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A46").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


		EwakstarDoer.State PL = new EwakstarDoer.State(
				"PL",
				new SimpleConfSpace.Builder()
						.addStrand(protein)
						.addStrand(ligand)
						.build()
		);

		int orderMag = 20; //order of magnitude worse in partition function we want to keep sequences relative to the wild-type
		double ewakstarEw = 20.0; //energy window within the wild-type for finding sequences in the "sequence filter" portion of ewakstarDoer
		int numPfConfs = 500; //num of conformations for the partition function calculation
		double pfEw = 1.0; //partition function energy window calculation
		int numTopOverallSeqs = 5; //end result number of sequences we want K* estimates for
		int numCpus = 4;
		double epsilon = 0.99;

		EwakstarDoer ewakstarDoer = new EwakstarDoer.Builder()
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
		PL.confTreeFactoryMin =(rcs) -> new ConfAStarTree.Builder(PL.emat, rcs)
				.setMaxNumNodes(2000000)
				.setTraditional()
				.build();

		PL.confTreeFactoryRigid =(rcs) -> new ConfAStarTree.Builder(PL.ematRigid, rcs)
				.setMaxNumNodes(2000000)
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
		double ewakstarEw = 10.0; //energy window within the wild-type for finding sequences in the "sequence filter" portion of ewakstarDoer
		int numPfConfs = 500; //num of conformations for the partition function calculation
		double pfEw = 1.0; //partition function energy window calculation
		int numTopOverallSeqs = 5; //end result number of sequences we want K* estimates for
		int numCpus = 4;
		double epsilon = 0.68;

		EwakstarDoer ewakstarDoer = new EwakstarDoer.Builder()
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
		PL.confTreeFactoryMin =(rcs) -> new ConfAStarTree.Builder(PL.emat, rcs)
				.setMaxNumNodes(2000000)
				.setTraditional()
				.build();

		PL.confTreeFactoryRigid =(rcs) -> new ConfAStarTree.Builder(PL.ematRigid, rcs)
				.setMaxNumNodes(2000000)
				.setTraditional()
				.build();

		return ewakstarDoer;
	}

}
