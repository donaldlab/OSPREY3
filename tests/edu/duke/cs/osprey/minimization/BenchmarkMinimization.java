package edu.duke.cs.osprey.minimization;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.gpu.GpuQueuePool;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;

public class BenchmarkMinimization extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
		// make a search problem
		System.out.println("Building search problem...");
		
		//String aaNames = "ALA VAL LEU ILE PHE TYR TRP CYS MET SER THR LYS ARG HIE HID ASP GLU ASN GLN GLY";
		//String aaNames = "ALA VAL LEU ILE";
		String aaNames = "ALA";
		//String mutRes = "39";
		String mutRes = "39 43";
		//String flexRes = "";
		//String flexRes = "40";
		//String flexRes = "40 41";
		String flexRes = "40 41 42 44 45";
		ArrayList<String> flexResList = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (String res : mutRes.split(" ")) {
			if (!res.isEmpty()) {
				flexResList.add(res);
				allowedAAs.add(new ArrayList<>(Arrays.asList(aaNames.split(" "))));
			}
		}
		for (String res : flexRes.split(" ")) {
			if (!res.isEmpty()) {
				flexResList.add(res);
				allowedAAs.add(new ArrayList<>());
			}
		}
		boolean doMinimize = true;
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = false;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"test", "test/1CC8/1CC8.ss.pdb", 
			flexResList, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion,
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null
		);
		
		final boolean useGpu = true;
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		if (useGpu) {
			egen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuQueuePool(1, 1, false));
		}
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		
		int numConfs = 10;
		
		// get the energy matrix
		File ematFile = new File("/tmp/emat.benchmarkMinimization.dat");
		if (ematFile.exists()) {
			System.out.println("\nReading energy matrix...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
		}
		if (search.emat == null) {
			System.out.println("\nComputing energy matrix...");
			search.emat =ecalc.calcEnergyMatrix(); 
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// get a few arbitrary conformations
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		List<ScoredConf> confs = new ArrayList<>();
		for (int i=0; i<numConfs; i++) {
			confs.add(tree.nextConf());
		}
		
		// build the energy function
		System.out.println("\nBuilding energy function...");
		EnergyFunction efunc = egen.fullConfEnergy(search.confSpace, search.shellResidues);
		
		// what do we expect the energies to be?
		double[] expectedEnergies = {
			-107.0147146474797,
			-107.14427781822944,
			-106.79145713460538,
			-106.92139365966113,
			-106.28769309286602,
			-107.11801397196132,
			-106.6789221601205,
			-106.41908247287255,
			-106.89600279647362,
			-106.74468003272366
		};
		
		// benchmark minimization
		System.out.println("\nBenchmarking...");
		Stopwatch cpuStopwatch = new Stopwatch().start();
		for (int i=0; i<confs.size(); i++) {
			
			ScoredConf conf = confs.get(i);
			
			RCTuple tuple = new RCTuple(conf.getAssignments());
			MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, search.confSpace, tuple);
	
			double energy = 0;
			
			final int simulateMinimization = 0;
			if (simulateMinimization > 0) {
				
				// don't actually minimize, just move the protein and calculate energies
				int numDofs = mof.getNumDOFs(); // is 17
				DoubleMatrix1D x = DoubleFactory1D.dense.make(numDofs);
				for (int d=0; d<numDofs; d++) {
					x.set(d, (mof.getConstraints()[0].get(d) + mof.getConstraints()[1].get(d))/2);
				}
				mof.setDOFs(x);
				
				if (simulateMinimization == 1) {
					
					// just slam the whole energy function a bunch
					for (int j=0; j<160; j++) {
						energy = mof.getValue(x);
					}
					
				} else if (simulateMinimization == 2) {
					
					// emulate CCD, move just one dof at a time
					for (int iter=0; iter<10; iter++) {
						for (int d=0; d<numDofs; d++) {
							for (int j=0; j<6; j++) {
								energy = mof.getValForDOF(d, x.get(d));
							}
						}
					}
				}
			
			} else {
				
				// minimize!
				Minimizer minimizer = new CCDMinimizer(mof, true);
				//Minimizer minimizer = new SimpleCCDMinimizer(mof);
				minimizer.minimize();
				energy = efunc.getEnergy();
			}
			
			// check the energy
			final double Epsilon = 1e-6;
			double absErr = energy - expectedEnergies[i];
			double relErr = absErr/Math.abs(expectedEnergies[i]);
			System.out.print(String.format("\texp:%12.8f  obs:%12.8f  absErr:%12.8f  relErr:%12.8f",
				expectedEnergies[i], energy,
				absErr, relErr
			));
			if (relErr > Epsilon) {
				System.out.print("     --=={  ENERGY TOO HIGH  }==--");
			}
			System.out.println();
		}
		cpuStopwatch.stop();
		System.out.println("finished in " + cpuStopwatch.getTime(2));
		
		/* what's the fastest speedup we can expect by optimizing the minimizer?
		
			the GPU energy function is about 15x faster than the CPU energy function on these energy functions
			but the minimizer only sees about a 2.5x speedup with the GPU energy function
			which means there's a lot of overhead in the minimizer itself!
			this benchmark takes the CCD minimizer about 7 seconds to finish with the CPU energy function
			and about 2.7 seconds with the GPU energy function
			which means if we properly optimize the minimizer itself to catch up to the energy function evaluation speed,
			we should be targeting about 0.5 s total run time for this benchmark with the GPU energy function
		*/
	}
}
