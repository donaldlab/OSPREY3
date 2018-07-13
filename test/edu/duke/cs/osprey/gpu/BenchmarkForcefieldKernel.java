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

package edu.duke.cs.osprey.gpu;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.opencl.kernels.ForcefieldKernelOpenCL;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.parallelism.TimingThread;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;
import edu.duke.cs.osprey.tupexp.LUTESettings;

@SuppressWarnings("unused")
public class BenchmarkForcefieldKernel extends TestBase {
	
	// NOTE: useful info for optimizing gpu kernels:
	// http://www.nvidia.com/content/GTC/documents/1068_GTC09.pdf
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
		// make a search problem
		System.out.println("Building search problem...");
		
		ResidueFlexibility resFlex = new ResidueFlexibility();
		resFlex.addMutable("39 43 46 47", "ALA");
		resFlex.addFlexible("40 41 42 44 45 48 49 50 51 52 53");
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
			"test", "examples/1CC8/1CC8.ss.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
			false, new ArrayList<>()
		);
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		//GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuQueuePool(1, 1, true));
		GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuStreamPool(1));
		
		benchmarkEfunc(search, egen, gpuegen);
		//benchmarkEmat(search, egen, gpuegen);
		//benchmarkMinimize(search, egen, gpuegen);
		//benchmarkEfuncFine(search.confSpace.m, egen, gpuegen);
	}
	
	private static void benchmarkEfunc(SearchProblem search, EnergyFunctionGenerator egen, GpuEnergyFunctionGenerator gpuegen)
	throws Exception {
		
		// figure out how many threads to use based on the gpu queue pool
		List<Integer> numThreadsList = new ArrayList<>();
		if (gpuegen.getOpenclQueuePool() != null) {
			for (int i=1; i<=gpuegen.getOpenclQueuePool().getNumQueues(); i<<=1) {
				numThreadsList.add(i);
			}
		} else {
			numThreadsList.add(1);
		}
		
		System.out.println("\nFull conf energy:");
		benchmarkEfunc(
			1000,
		//profileEfunc(
		//	5000,
			search.confSpace.m,
			new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return egen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
				}
			},
			new Factory<GpuForcefieldEnergy,Molecule>() {
				@Override
				public GpuForcefieldEnergy make(Molecule mol) {
					return gpuegen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
				}
			},
			numThreadsList
		);
		
		System.out.println("\nIntra and shell energy:");
		benchmarkEfunc(
		//profileEfunc(
			40000,
			search.confSpace.m,
			new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					Residue res = search.confSpace.posFlex.get(0).res;
					return egen.intraAndShellEnergy(getResidue(res, mol), getResidues(search.shellResidues, mol)); 
				}
			},
			new Factory<GpuForcefieldEnergy,Molecule>() {
				@Override
				public GpuForcefieldEnergy make(Molecule mol) {
					Residue res = search.confSpace.posFlex.get(0).res;
					return gpuegen.intraAndShellEnergy(getResidue(res, mol), getResidues(search.shellResidues, mol));
				}
			},
			numThreadsList
		);
		
		System.out.println("\nPairwise energy:");
		// GPU is actually significantly slower for these terms
		// transfers between the CPU and GPU cause a LOT of overhead!
		// I've hammered on this for a while, but I haven't found any faster ways of doing the transfers
		// mapping buffers and pinning memory seems to be significantly slower than direct transfers
		// at least on my hardware at home =(
		Residue res1 = search.confSpace.posFlex.get(0).res; // GLY
		Residue res2 = search.confSpace.posFlex.get(2).res; // GLY
		//Residue res1 = search.confSpace.m.getResByPDBResNumber("65"); // LYS
		//Residue res2 = search.confSpace.m.getResByPDBResNumber("68"); // ARG
		benchmarkEfunc(
		//profileEfunc(
			50000,
			search.confSpace.m,
			new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return egen.resPairEnergy(getResidue(res1, mol), getResidue(res2, mol));
				}
			},
			new Factory<GpuForcefieldEnergy,Molecule>() {
				@Override
				public GpuForcefieldEnergy make(Molecule mol) {
					return gpuegen.resPairEnergy(getResidue(res1, mol), getResidue(res2, mol));
				}
			},
			numThreadsList
		);
	}
	
	private static Residue getResidue(Residue res, Molecule mol) {
		return mol.residues.get(res.indexInMolecule);
	}
	
	private static List<Residue> getResidues(List<Residue> residues, Molecule mol) {
		List<Residue> matched = new ArrayList<>();
		for (Residue res : residues) {
			matched.add(getResidue(res, mol));
		}
		return matched;
	}
	
	private static void benchmarkEfunc(int numRuns, Molecule baseMol, Factory<EnergyFunction,Molecule> efuncs, Factory<GpuForcefieldEnergy,Molecule> gpuefuncs, List<Integer> numThreadsList) {
		
		// get the answer
		final double expectedEnergy = efuncs.make(baseMol).getEnergy();
		
		// benchmark cpus
		Stopwatch cpuStopwatch = null;
		for (int numThreads : numThreadsList) {
			System.out.print("Benchmarking " + numThreads + " CPUs... ");
			
			// make the thread pool
			List<TimingThread> threads = new ArrayList<>();
			for (int i=0; i<numThreads; i++) {
				threads.add(new TimingThread("Cpu-" + i) {
					
					private Molecule mol;
					private EnergyFunction efunc;
					
					@Override
					public void warmup() {
						mol = new Molecule(baseMol);
						efunc = efuncs.make(mol);
						for (int k=0; k<10; k++) {
							efunc.getEnergy();
						}
					}
					
					@Override
					public void time() {
						
						double energy = 0;
						int numLocalRuns = numRuns/numThreads;
						for (int k=0; k<numLocalRuns; k++) {
							energy = efunc.getEnergy();
						}
						
						checkEnergy(expectedEnergy, energy);
					}
				});
			}
			
			// run threads and wait
			Stopwatch cpusStopwatch = TimingThread.timeThreads(threads);
			
			if (numThreads == 1) {
				cpuStopwatch = cpusStopwatch;
			}
			
			// show timing info
			System.out.println(String.format(" finished in %s, avg time per op: %s, speedup: %.2fx",
				cpusStopwatch.getTime(2),
				TimeFormatter.format(cpusStopwatch.getTimeNs()/numRuns, TimeUnit.MICROSECONDS),
				(double)cpuStopwatch.getTimeNs()/cpusStopwatch.getTimeNs()
			));
		}
		
		// benchmark gpus
		for (int numThreads : numThreadsList) {
			System.out.print("Benchmarking " + numThreads + " GPUs... ");
			
			// make the thread pool
			List<TimingThread> threads = new ArrayList<>();
			for (int i=0; i<numThreads; i++) {
				threads.add(new TimingThread("Gpu-" + i) {
					
					private Molecule mol;
					private GpuForcefieldEnergy gpuefunc;
					
					@Override
					public void warmup() {
						mol = new Molecule(baseMol);
						gpuefunc = gpuefuncs.make(mol);
						for (int j=0; j<10; j++) {
							gpuefunc.getEnergy();
						}
					}
					
					@Override
					public void time() {
						try {
							
							double energy = 0;
							int numLocalRuns = numRuns/numThreads;
							for (int k=0; k<numLocalRuns; k++) {
								energy = gpuefunc.getEnergy();
							}
							
							checkEnergy(expectedEnergy, energy);
						
						} finally {
							gpuefunc.clean();
						}
					}
				});
			}
			
			// run the threads and wait
			Stopwatch gpuStopwatch = TimingThread.timeThreads(threads);
			
			// show timing info
			System.out.println(String.format(" finished in %s, avg time per op: %s, speedup: %.2fx",
				gpuStopwatch.getTime(2),
				TimeFormatter.format(gpuStopwatch.getTimeNs()/numRuns, TimeUnit.MICROSECONDS),
				(double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs()
			));
		}
	}
	
	private static void profileEfunc(int numRuns, Molecule baseMol, Factory<EnergyFunction,Molecule> efuncs, Factory<GpuForcefieldEnergy,Molecule> gpuefuncs, List<Integer> numThreadsList) {
		
		GpuForcefieldEnergy gpuefunc = gpuefuncs.make(baseMol);
		
		// get the answer
		final double expectedEnergy = efuncs.make(baseMol).getEnergy();
		
		// do warmup
		for (int i=0; i<numRuns/10; i++) {
			gpuefunc.getEnergy();
		}
		
		ForcefieldKernelOpenCL kernel = (ForcefieldKernelOpenCL)gpuefunc.getKernel();
		
		Stopwatch stopwatch = new Stopwatch().start();
		
		// do profiling
		double energy = 0;
		for (int i=0; i<numRuns; i++) {
			energy = gpuefunc.getEnergy();
		}
		
		stopwatch.stop();
		
		// print the report
		System.out.println("GPU profiling info:");
		System.out.println("atom pairs:      " + kernel.getForcefield().getFullSubset().getNumAtomPairs());
		System.out.println("GPU memory used: " + kernel.getGpuBytesNeeded()/1024 + " KiB");
		System.out.println("us per op: " + TimeFormatter.format(stopwatch.getTimeNs()/numRuns, TimeUnit.MICROSECONDS));
		
		checkEnergy(expectedEnergy, energy);
		
		gpuefunc.clean();
	}
	
	protected static void checkEnergy(double exp, double obs) {
		final double Epsilon = 1e-12;
		double absErr = Math.abs(exp - obs);
		double relErr = absErr/Math.abs(exp);
		if (relErr > Epsilon) {
			System.out.println(String.format("Wrong energy! exp: %12.6f  obs: %12.6f  absErr: %12.6f  relErr: %12.6f",
				exp, obs, absErr, relErr
			));
		}
	}

	private static void benchmarkMinimize(SearchProblem search, EnergyFunctionGenerator egen, GpuEnergyFunctionGenerator gpuegen)
	throws Exception {
		
		int numConfs = 16;
		
		// get a few arbitrary conformations
		search.emat = new SimpleEnergyMatrixCalculator.Cpu(2, egen.ffParams, search.confSpace, search.shellResidues).calcEnergyMatrix();
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setMPLP(new ConfAStarTree.MPLPBuilder()
				.setNumIterations(4)
			).build();
		List<ScoredConf> confs = new ArrayList<>();
		for (int i=0; i<numConfs; i++) {
			confs.add(tree.nextConf());
		}
		
		double energy;
		
		// benchmark the cpu
		System.out.println("\nBenchmarking CPU...");
		Stopwatch cpuStopwatch = new Stopwatch().start();
		EnergyFunction efunc = egen.fullConfEnergy(search.confSpace, search.shellResidues);
		energy = 0;
		for (ScoredConf conf : confs) {
			energy += minimize(efunc, search.confSpace, conf);
		}
		cpuStopwatch.stop();
		System.out.println("\tfinished in " + cpuStopwatch.getTime(2));
		System.out.println("\tenergy sum: " + energy);
		
		// benchmark the gpu
		System.out.println("\nBenchmarking GPU...");
		Stopwatch gpuStopwatch = new Stopwatch().start();
		GpuForcefieldEnergy gpuefunc = gpuegen.fullConfEnergy(search.confSpace, search.shellResidues);
		energy = 0;
		for (ScoredConf conf : confs) {
			energy += minimize(gpuefunc, search.confSpace, conf);
		}
		gpuStopwatch.stop();
		gpuefunc.clean();
		System.out.println("\tfinished in " + gpuStopwatch.getTime(2));
		System.out.println("\tenergy sum: " + energy);
		
		// calculate speedup
		System.out.println(String.format("\nspeedup: %.2fx", (double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs())); 
	}
	
	private static double minimize(EnergyFunction efunc, ConfSpace confSpace, ScoredConf conf) {
		RCTuple tuple = new RCTuple(conf.getAssignments());
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple);
		new CCDMinimizer(mof, true).minimize();
		//new SimpleCCDMinimizer(mof).minimize();
		return efunc.getEnergy();
	}
	
	private static void benchmarkEfuncFine(Molecule mol, EnergyFunctionGenerator egen, GpuEnergyFunctionGenerator gpuegen) {
		
		List<Integer> threadsList = Arrays.asList(1);
		
		class Run {
			
			public int numResidues;
			public int numRuns;
			
			public Run(int numResidues, int numRuns) {
				this.numResidues = numResidues;
				this.numRuns = numRuns;
			}
		}
		
		Deque<Run> runs = new ArrayDeque<>();
		runs.add(new Run(1, 1000000));
		runs.add(new Run(2, 500000));
		runs.add(new Run(3, 200000));
		runs.add(new Run(4, 100000));
		runs.add(new Run(5, 80000));
		runs.add(new Run(6, 60000));
		runs.add(new Run(7, 40000));
		runs.add(new Run(8, 30000));
		runs.add(new Run(9, 20000));
		runs.add(new Run(10, 16000));
		runs.add(new Run(11, 14000));
		runs.add(new Run(12, 13000));
		runs.add(new Run(13, 12000));
		runs.add(new Run(14, 11500));
		runs.add(new Run(15, 11000));
		runs.add(new Run(16, 10500));
		runs.add(new Run(17, 10000));
		runs.add(new Run(18, 9500));
		runs.add(new Run(19, 9000));
		runs.add(new Run(20, 8500));
		runs.add(new Run(25, 7000));
		runs.add(new Run(30, 6000));
		runs.add(new Run(35, 4000));
		runs.add(new Run(40, 2000));
		runs.add(new Run(45, 1000));
		runs.add(new Run(50, 900));
		runs.add(new Run(55, 800));
		runs.add(new Run(60, 700));
		runs.add(new Run(65, 650));
		runs.add(new Run(70, 620));
		
		//final Run run = runs.getLast();
		
		for (Run run : runs) {
			for (int i=0; i<3; i++) {
				benchmarkEfunc(
					run.numRuns,
					mol,
					new Factory<EnergyFunction,Molecule>() {
						@Override
						public EnergyFunction make(Molecule mol) {
							return makeFineEfunc(mol, egen, run.numResidues);
						}
					},
					new Factory<GpuForcefieldEnergy,Molecule>() {
						@Override
						public GpuForcefieldEnergy make(Molecule mol) {
							return makeFineGpuEfunc(mol, gpuegen, run.numResidues);
						}
					},
					threadsList
				);
			}
		}
	}
	
	private static EnergyFunction makeFineEfunc(Molecule mol, EnergyFunctionGenerator egen, int numResidues) {
		
		MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
		
		for (int pos1=0; pos1<numResidues; pos1++) {
			
			Residue res1 = mol.residues.get(pos1);
			efunc.addTerm(new SingleResEnergy(res1, egen.ffParams));
			
			for (int pos2=0; pos2<pos1; pos2++) {
				
				Residue res2 = mol.residues.get(pos2);
				efunc.addTerm(new ResPairEnergy(res1, res2, egen.ffParams));
			}
		}
		
		return efunc;
	}
	
	private static GpuForcefieldEnergy makeFineGpuEfunc(Molecule mol, GpuEnergyFunctionGenerator egen, int numResidues) {
		
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		
		for (int pos1=0; pos1<numResidues; pos1++) {
			
			Residue res1 = mol.residues.get(pos1);
			interactions.addResidue(res1);
			
			for (int pos2=0; pos2<pos1; pos2++) {
				
				Residue res2 = mol.residues.get(pos2);
				interactions.addResiduePair(res1, res2);
			}
		}
		
		return new GpuForcefieldEnergy(egen.ffParams, interactions, egen.getOpenclQueuePool());
	}
}
