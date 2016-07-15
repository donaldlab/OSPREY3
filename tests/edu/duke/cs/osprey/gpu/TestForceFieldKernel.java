package edu.duke.cs.osprey.gpu;

import java.util.List;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.kernels.ForceFieldKernel;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
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
		// cheat sheet:  15 GLY, 17 SER, 18 TRP, 25 TRP, 22 ARG, 24 ALA
		
		// small pair GLY-SER
		//List<Atom> atomsA = mol.getResByPDBResNumber("15").atoms;
		//List<Atom> atomsB = mol.getResByPDBResNumber("17").atoms;
		
		// large pair TRP-TRP
		List<Atom> atomsA = mol.getResByPDBResNumber("18").atoms;
		List<Atom> atomsB = mol.getResByPDBResNumber("25").atoms;
		
		System.out.println("num atom pairs: " + atomsA.size()*atomsB.size());
		
		// make the energy function
		ForcefieldParams ffparams = TestBase.makeDefaultFFParams();
		ForcefieldEnergy ffenergy = new ForcefieldEnergy(false, atomsA, atomsB, ffparams);
		ffenergy.setUseCache(false);
		int numTerms = ffenergy.getNumTerms();
		System.out.println("energy terms: " + numTerms);
		
		final int NumRuns = 10000;
		
		// benchmark the cpu
		double energy = 0;
		System.out.println("running forcefield on CPU...");
		Stopwatch cpuStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			energy = ffenergy.calculateTotalEnergy();
		}
		cpuStopwatch.stop();
		System.out.println("Cpu time: " + cpuStopwatch.getTime(1));
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
		System.out.println("running forcefield on GPU...");
		Stopwatch gpuStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			kernel.uploadRunDownloadSync();
			//kernel.runAsync();
		}
		kernel.waitForGpu();
		energy = kernel.getOut().get(0);
		gpuStopwatch.stop();
		System.out.println("Gpu time: " + gpuStopwatch.getTime(1));
		System.out.println("Energy: " + energy);
		
		System.out.println(String.format("%.2fx speedup", (float)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs()));
	}
	
	private static void checkTotalEnergy(Molecule mol) {
		
		System.out.println("Building energy function...");
		EnergyFunction efunc = EnvironmentVars.curEFcnGenerator.fullMolecEnergy(mol);
		
		final int NumRuns = 10;
		
		// benchmark the cpu
		double energy = 0;
		System.out.println("running forcefield on CPU...");
		Stopwatch cpuStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			energy = efunc.getEnergy();
		}
		System.out.println("Cpu time: " + cpuStopwatch.getTime(1));
		System.out.println("Energy: " + energy);
		
		// TODO: compare to gpu
	}
}
