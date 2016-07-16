package edu.duke.cs.osprey.gpu;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.gpu.kernels.ForceFieldKernel;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Stopwatch;

public class TestForceFieldKernel extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		System.out.println("Initializing...");
		initDefaultEnvironment();
		MultiTermEnergyFunction.setNumThreads(1); // turn off multi-threading
		
		// get a molecule
		System.out.println("Reading PDB...");
		Molecule mol = PDBFileReader.readPDBFile("test/DAGK/2KDC.P.forOsprey.pdb");
		
		// pick some subset of the molecule for energy calculations
		Residue gly06 = mol.getResByPDBResNumber("6");
		Residue gly15 = mol.getResByPDBResNumber("15");
		Residue ser17 = mol.getResByPDBResNumber("17");
		Residue trp18 = mol.getResByPDBResNumber("18");
		Residue trp25 = mol.getResByPDBResNumber("25");
		Residue arg22 = mol.getResByPDBResNumber("22");
		Residue ala24 = mol.getResByPDBResNumber("24");
		Residue ile26 = mol.getResByPDBResNumber("26");
		Residue phe31 = mol.getResByPDBResNumber("31");
		Residue arg32 = mol.getResByPDBResNumber("32");
		Residue glu34 = mol.getResByPDBResNumber("34");
		Residue val36 = mol.getResByPDBResNumber("36");
		Residue leu39 = mol.getResByPDBResNumber("39");
		Residue trp47 = mol.getResByPDBResNumber("47");
		Residue leu48 = mol.getResByPDBResNumber("48");
		Residue ile53 = mol.getResByPDBResNumber("53");
		Residue arg55 = mol.getResByPDBResNumber("55");
		Residue val56 = mol.getResByPDBResNumber("56");
		Residue leu57 = mol.getResByPDBResNumber("57");
		Residue ile59 = mol.getResByPDBResNumber("59");
		Residue val62 = mol.getResByPDBResNumber("62");
		Residue leu64 = mol.getResByPDBResNumber("64");
		Residue val65 = mol.getResByPDBResNumber("65");
		Residue met66 = mol.getResByPDBResNumber("66");
		
		//Residue[] residues = { gly15 };
		//Residue[] residues = { gly06, gly15 };
		//Residue[] residues = { gly15, ser17 }; // 2x speedup
		//Residue[] residues = { gly15, ser17, trp18, trp25 }; // 19x speedup
		//Residue[] residues = { gly15, ser17, trp18, trp25, arg22, ala24 }; // 41x speedup
		//Residue[] residues = { gly15, ser17, trp18, trp25, arg22, ala24, ile26, phe31, arg32, glu34 }; // 80x speedup
		//Residue[] residues = { gly15, ser17, trp18, trp25, arg22, ala24, ile26, phe31, arg32, glu34, val36, leu39, trp47, leu48 }; // 95x speedup
		Residue[] residues = { gly15, ser17, trp18, trp25, arg22, ala24, ile26, phe31, arg32, glu34, val36, leu39, trp47, leu48,
			ile53, arg55, val56, leu57, ile59, val62, leu64, val65, met66 }; // 92x speedup
		
		System.out.println("\nBuilding energy function...");
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		int numTerms = 0;
		for (int pos1=0; pos1<residues.length; pos1++) {
			
			SingleResEnergy single = new SingleResEnergy(residues[pos1], egen.ffParams);
			numTerms += single.getFFEnergy().getNumTerms();
			efunc.addTerm(single);
			interactions.addResidue(residues[pos1]);
			
			for (int pos2=0; pos2<pos1; pos2++) {
				
				ResPairEnergy pair = new ResPairEnergy(residues[pos1], residues[pos2], egen.ffParams);
				numTerms += pair.getFFEnergy().getNumTerms();
				efunc.addTerm(pair);
				interactions.addResiduePair(residues[pos1], residues[pos2]);
			}
		}
		System.out.println("energy terms: " + numTerms);
		
		// make the big forcefield energy function
		System.out.println("\nBuilding big forcefield...");
		
		ForcefieldParams ffparams = TestBase.makeDefaultFFParams();
		BigForcefieldEnergy ffenergy = new BigForcefieldEnergy(ffparams, interactions);
		System.out.println("energy terms: " + ffenergy.getNumAtomPairs());
		
		final int NumRuns = 1000;
		
		// benchmark the cpu on the energy function
		double energy = 0;
		System.out.println("\nrunning energy function on CPU...");
		Stopwatch cpuEfuncStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			energy = efunc.getEnergy();
		}
		cpuEfuncStopwatch.stop();
		System.out.println("Cpu time: " + cpuEfuncStopwatch.getTime(2));
		System.out.println("Energy: " + energy);
		
		// benchmark the cpu on the big forcefield
		System.out.println("\nrunning big forcefield on CPU...");
		Stopwatch cpuBigffStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			energy = ffenergy.calculateTotalEnergy();
		}
		cpuBigffStopwatch.stop();
		System.out.println("Cpu time: " + cpuBigffStopwatch.getTime(2));
		System.out.println("Energy: " + energy);
		
		// prep the gpu
		ForceFieldKernel.Bound kernel = new ForceFieldKernel().bind();
		kernel.setWorkSize(numTerms);
		for (int i=0; i<numTerms; i++) {
			kernel.getIn().put(i);
		}
		kernel.getIn().rewind();
		kernel.uploadSync();
		
		// benchmark the gpu
		System.out.println("\nrunning forcefield on GPU...");
		Stopwatch gpuStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			//kernel.uploadRunDownloadSync();
			kernel.runAsync();
		}
		kernel.waitForGpu();
		energy = kernel.getOut().get(0);
		gpuStopwatch.stop();
		System.out.println("Gpu time: " + gpuStopwatch.getTime(2));
		System.out.println("Energy: " + energy);
		
		System.out.println(String.format("%.2fx speedup", (float)cpuBigffStopwatch.getTimeNs()/gpuStopwatch.getTimeNs()));
	}
}
