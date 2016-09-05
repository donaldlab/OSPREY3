package edu.duke.cs.osprey.gpu;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
import edu.duke.cs.osprey.confspace.ConfSpace;
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
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class BenchmarkForcefieldKernel extends TestBase {
	
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
			flexResList, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null, false
		);
		
		//benchmarkEfunc(search);
		//benchmarkEmat(search);
		benchmarkMinimize(search);
	}
	
	private static void benchmarkEfunc(SearchProblem search)
	throws Exception {
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams());
		
		System.out.println("NOTE: disable the energy cache in ForcefieldEnergy, or these tests will make the GPU looks really bad! =P");
		
		System.out.println("\nFull conf energy:");
		benchmarkEfunc(1000,
			egen.fullConfEnergy(search.confSpace, search.shellResidues),
			gpuegen.fullConfEnergy(search.confSpace, search.shellResidues)
		);
		
		System.out.println("\nIntra and shell energy:");
		benchmarkEfunc(6000,
			egen.intraAndShellEnergy(search.confSpace.posFlex.get(0).res, search.shellResidues),
			gpuegen.intraAndShellEnergy(search.confSpace.posFlex.get(0).res, search.shellResidues)
		);
		
		System.out.println("\nPairwise energy:");
		// TODO: GPU is actually significantly slower for these terms
		// there's so few atom pairs, the overhead with the GPU is slowing us down
		// need to optimize more, maybe look into faster memory transfers?
		benchmarkEfunc(100000,
			egen.resPairEnergy(search.confSpace.posFlex.get(0).res, search.confSpace.posFlex.get(2).res),
			gpuegen.resPairEnergy(search.confSpace.posFlex.get(0).res, search.confSpace.posFlex.get(2).res)
		);
	}
	
	private static void benchmarkEfunc(int numRuns, EnergyFunction efunc, GpuForcefieldEnergy gpuefunc) {
		
		// benchmark the cpu
		System.out.print("Benchmarking CPU...");
		Stopwatch cpuStopwatch = new Stopwatch().start(); 
		for (int i=0; i<numRuns; i++) {
			efunc.getEnergy();
		}
		cpuStopwatch.stop();
		System.out.println(" finished in " + cpuStopwatch.getTime(2));
		
		// benchmark the gpu
		System.out.print("Benchmarking GPU...");
		Stopwatch gpuStopwatch = new Stopwatch().start();
		for (int i=0; i<numRuns; i++) {
			gpuefunc.getEnergy();
		}
		gpuStopwatch.stop();
		gpuefunc.cleanup();
		System.out.println(String.format(" finished in %s, speedup: %.2fx, numPairs: %d, GPU mem used: %.2f MiB",
			gpuStopwatch.getTime(2),
			(double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs(),
			gpuefunc.getForcefieldEnergy().getNumAtomPairs(),
			(double)gpuefunc.getGpuBytesNeeded()/1024/1024
		)); 
	}
	
	private static void benchmarkEmat(SearchProblem search)
	throws Exception {
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		
		GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams());
		SimpleEnergyCalculator gpuecalc = new SimpleEnergyCalculator(gpuegen, search.confSpace, search.shellResidues);
		
		// benchmark the cpu
		System.out.println("\nBenchmarking CPU...");
		Stopwatch cpuStopwatch = new Stopwatch().start();
		EnergyMatrix emat = ecalc.calcEnergyMatrix();
		cpuStopwatch.stop();
		
		// benchmark the gpu
		System.out.println("\nBenchmarking GPU...");
		Stopwatch gpuStopwatch = new Stopwatch().start();
		EnergyMatrix gpuemat = gpuecalc.calcEnergyMatrix();
		gpuStopwatch.stop();
		
		// calculate speedup
		System.out.println(String.format("\nspeedup: %.2fx", (double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs()));
		
		// check the result
		// TODO: apparently these results don't match, need to find out why
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
				
				assertThat(emat.getOneBody(pos1, rc1), isRelatively(gpuemat.getOneBody(pos1, rc1)));
				
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++) {
						
						assertThat(emat.getPairwise(pos1, rc1, pos2, rc2), isRelatively(gpuemat.getPairwise(pos1, rc1, pos2, rc2)));
					}
				}
			}
		}
	}
	
	private static void benchmarkMinimize(SearchProblem search)
	throws Exception {
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		
		GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams());
		
		int numConfs = 10;
		
		// get a few arbitrary conformations
		search.emat = ecalc.calcEnergyMatrix();
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		List<ScoredConf> confs = new ArrayList<>();
		for (int i=0; i<numConfs; i++) {
			confs.add(tree.nextConf());
		}
		
		double energy;
		
		// benchmark the cpu
		System.out.println("\nBenchmarking CPU...");
		EnergyFunction efunc = egen.fullConfEnergy(search.confSpace, search.shellResidues);
		Stopwatch cpuStopwatch = new Stopwatch().start();
		energy = 0;
		for (ScoredConf conf : confs) {
			energy += minimize(efunc, search.confSpace, conf);
		}
		cpuStopwatch.stop();
		System.out.println("\tfinished in " + cpuStopwatch.getTime(2));
		System.out.println("\tenergy sum: " + energy);
		
		// benchmark the gpu
		System.out.println("\nBenchmarking GPU...");
		GpuForcefieldEnergy gpuefunc = gpuegen.fullConfEnergy(search.confSpace, search.shellResidues);
		Stopwatch gpuStopwatch = new Stopwatch().start();
		energy = 0;
		for (ScoredConf conf : confs) {
			energy += minimize(gpuefunc, search.confSpace, conf);
		}
		gpuStopwatch.stop();
		gpuefunc.cleanup();
		System.out.println("\tfinished in " + gpuStopwatch.getTime(2));
		System.out.println("\tenergy: " + energy);
		
		// calculate speedup
		System.out.println(String.format("\nspeedup: %.2fx", (double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs())); 
	}
	
	private static double minimize(EnergyFunction efunc, ConfSpace confSpace, ScoredConf conf) {
		RCTuple tuple = new RCTuple(conf.getAssignments());
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple);
		new CCDMinimizer(mof, true).minimize();
		return efunc.getEnergy();
	}
}
