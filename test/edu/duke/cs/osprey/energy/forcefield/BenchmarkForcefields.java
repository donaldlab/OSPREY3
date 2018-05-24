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

package edu.duke.cs.osprey.energy.forcefield;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.ResidueForcefieldEnergyCuda;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;

public class BenchmarkForcefields extends TestBase {
	
	private static class Result {
		
		double ms1;
		double ms2;
		double ms3;
		
		public Result(double ms1, double ms2, double ms3) {
			this.ms1 = ms1;
			this.ms2 = ms2;
			this.ms3 = ms3;
		}
	}
	
	public static void main(String[] args) {
		
		ForcefieldParams ffparams = new ForcefieldParams();
		
		// get a conf space
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
		strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
		// pre-compute atom connectivities
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(confSpace)
			.setParallelism(Parallelism.makeCpu(4))
			.build();
		ResPairCache resPairCache = new ResPairCache(ffparams, connectivity);
		
		// get a molecule
		Molecule mol = confSpace.makeMolecule(new int[] { 0, 0, 0, 0, 0 }).mol;
		
		benchmarkComparison(mol, confSpace, ffparams, resPairCache);
		//benchmarkGpu(mol, confSpace, ffparams, resPairCache);
	}
	
	private static void benchmarkComparison(Molecule mol, SimpleConfSpace confSpace, ForcefieldParams ffparams, ResPairCache resPairCache) {
		
		// multi term energy function
		Result base = benchmark(MultiTermEnergyFunction.class.getSimpleName(), null, () -> {
			MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
			for (int pos1=0; pos1<confSpace.positions.size(); pos1++) {
				Residue res1 = mol.getResByPDBResNumber(confSpace.positions.get(pos1).resNum);
				efunc.addTerm(new SingleResEnergy(res1, ffparams));
				for (int pos2=0; pos2<pos1; pos2++) {
					Residue res2 = mol.getResByPDBResNumber(confSpace.positions.get(pos2).resNum);
					efunc.addTerm(new ResPairEnergy(res1, res2, ffparams));
				}
				for (String shellResNum : confSpace.shellResNumbers) {
					Residue shellRes = mol.getResByPDBResNumber(shellResNum);
					efunc.addTerm(new ResPairEnergy(res1, shellRes, ffparams));
				}
			}
			return efunc;
		});
		
		RCTuple frag = new RCTuple(new int[] { 0, 0, 0, 0, 0 });
		ResidueInteractions inters = ResInterGen.of(confSpace)
			.addIntras(frag)
			.addInters(frag)
			.addShell(frag)
			.make();
		
		// residue forcefield
		benchmark(ResidueForcefieldEnergy.class.getSimpleName(), base, () -> {
			// new ResidueForcefieldEnergy(resPairCache, interactions, mol);
			return new ResidueForcefieldEnergy(resPairCache, inters, mol);
		});
		
		// residue forcefield cuda
		GpuStreamPool streams = new GpuStreamPool(1, 1);
		benchmark(ResidueForcefieldEnergyCuda.class.getSimpleName(), base, () -> {
			return new ResidueForcefieldEnergyCuda(streams, resPairCache, inters, mol);
		});
		streams.cleanup();
		
		// big forcefield
		benchmark(BigForcefieldEnergy.class.getSimpleName(), base, () -> {
			ForcefieldInteractions ffinters = new ForcefieldInteractions();
			for (int pos1=0; pos1<confSpace.positions.size(); pos1++) {
				Residue res1 = mol.getResByPDBResNumber(confSpace.positions.get(pos1).resNum);
				ffinters.addResidue(res1);
				for (int pos2=0; pos2<pos1; pos2++) {
					Residue res2 = mol.getResByPDBResNumber(confSpace.positions.get(pos2).resNum);
					ffinters.addResiduePair(res1, res2);
				}
				for (String shellResNum : confSpace.shellResNumbers) {
					Residue shellRes = mol.getResByPDBResNumber(shellResNum);
					ffinters.addResiduePair(res1, shellRes);
				}
			}
			return new BigForcefieldEnergy(ffparams, ffinters);
		});
	}
	
	private static Result benchmark(String name, Result base, Supplier<EnergyFunction> efuncs) {
		
		System.out.println(name);
		
		Result result = new Result(0, 0, 0);
		
		// create-cleanup cycles
		result.ms1 = benchmark("create-cleanup", 10, 200, 3, base == null ? null : base.ms1, () -> {
			EnergyFunction.Tools.cleanIfNeeded(efuncs.get());
		});
		
		// second benchmark, run cycles
		{
			EnergyFunction efunc = efuncs.get();
			try {
				result.ms2 = benchmark("run", 10, 600, 3, base == null ? null : base.ms2, () -> {
					efunc.getEnergy();
				});
			} finally {
				EnergyFunction.Tools.cleanIfNeeded(efunc);
			}
		}
		
		// first benchmark, create-run-cleanup cycles
		result.ms3 = benchmark("create-run-cleanup", 10, 100, 3, base == null ? null : base.ms3, () -> {
			EnergyFunction efunc = efuncs.get();
			try {
				efunc.getEnergy();
			} finally {
				EnergyFunction.Tools.cleanIfNeeded(efunc);
			}
		});
		
		
		return result;
	}
	
	private static double benchmark(String name, int numWarmupRuns, int numRuns, int replicates, Double baseMs, Runnable iter) {
		
		List<Double> ms = new ArrayList<>();
		
		for (int repl=0; repl<replicates; repl++) {
			
			// warm up first
			for (int i=0; i<numWarmupRuns; i++) {
				iter.run();
			}
			
			// time it for reals
			Stopwatch stopwatch = new Stopwatch().start();
			for (int i=0; i<numRuns; i++) {
				iter.run();
			}
			stopwatch.stop();
			
			// show the results
			System.out.print(String.format("\t%22s finished in %10s   %6.1f ops",
				name,
				stopwatch.stop().getTime(2),
				numRuns/stopwatch.getTimeS()
			));
			if (baseMs != null) {
				System.out.print(String.format("   %.2fx speedup", baseMs/stopwatch.getTimeMs()));
			}
			System.out.println();
			
			ms.add(stopwatch.getTimeMs());
		}
		
		// get the best results
		double bestMs = ms.stream()
			.reduce((Double a, Double b) -> Math.min(a, b))
			.get();
		
		System.out.print(String.format("\t                     best time: %10s   %6.1f ops",
			TimeFormatter.format((long)(bestMs*1000*1000), 2),
			(double)numRuns/bestMs*1000
		));
		if (baseMs != null) {
			System.out.print(String.format("   %.2fx speedup", baseMs/bestMs));
		}
		System.out.println();
		
		System.out.println();
		
		return bestMs;
	}
	
	private static void benchmarkGpu(Molecule mol, SimpleConfSpace confSpace, ResPairCache resPairCache) {
		
		RCTuple frag = new RCTuple(new int[] { 0, 0, 0, 0, 0 });
		ResidueInteractions inters = ResInterGen.of(confSpace)
			.addIntras(frag)
			.addInters(frag)
			.make();
		
		GpuStreamPool streams = new GpuStreamPool(1, 1);
		
		benchmark("create-cleanup", 100, 1000, 5, null, () -> {
			ResidueForcefieldEnergyCuda efunc = new ResidueForcefieldEnergyCuda(streams, resPairCache, inters, mol);
			efunc.getEnergy();
			efunc.clean();
		});
			
		/*
		try {
			ResidueForcefieldEnergyCuda efunc = new ResidueForcefieldEnergyCuda(streams, resPairCache, inters, mol);
			benchmark("run", 1000, 20000, 5, null, () -> {
				efunc.getEnergy();
			});
			efunc.cleanup();
		} catch (IOException ex) {
			throw new Error(ex);
		}
		*/
			
		streams.cleanup();
		
		// best so far for run: 7131.8 ops
	}
}
