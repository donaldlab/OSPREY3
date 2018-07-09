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

package edu.duke.cs.osprey;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.ParallelEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Stopwatch;

public class EnergyProfiling {
	
	@SuppressWarnings("unused")
	public static void main(String[] args) {
		
		// check the cwd
		String path = new File("").getAbsolutePath();
		if (!path.endsWith("examples/DAGK")) {
			throw new Error("This profiler was designed to run in the examples/DAGK folder\n\tcwd: " + path);
		}

		// load configuration
		ConfigFileParser cfp = new ConfigFileParser();
		cfp.loadData();
		
		// configure energy function parallelization
		final int NumThreads = 4;
		MultiTermEnergyFunction.setNumThreads(NumThreads);
		ParallelEnergyFunction.startCrew(NumThreads);
		
		// read a big test protein, the bigger the better
		Molecule mol = new Strand.Builder(PDBIO.readFile("2KDC.P.forOsprey.pdb")).build().mol;
		
		System.out.println("\n\nBuilding energy functions...");
		
		// get all the energy function terms
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		ArrayList<SingleResEnergy> singleTerms = new ArrayList<>();
		for (Residue res : mol.residues) {
			singleTerms.add(new SingleResEnergy(res, ffparams));
		}
		ArrayList<ResPairEnergy> pairTerms = new ArrayList<>();
		for (int i=0; i<mol.residues.size(); i++) {
			for (int j=0; j<i; j++) {
				pairTerms.add(new ResPairEnergy(mol.residues.get(i), mol.residues.get(j), ffparams));
			}
		}
		
		System.out.println(String.format("Built %d energy terms:\n\t%d singles\n\t%d pairs",
			singleTerms.size() + pairTerms.size(),
			singleTerms.size(),
			pairTerms.size()
		));

		// build a few different energy functions to test
		
		// the whole shebang
		MultiTermEnergyFunction totalEFunc = new MultiTermEnergyFunction();
		for (SingleResEnergy term : singleTerms) {
			totalEFunc.addTerm(term);
		}
		for (ResPairEnergy term : pairTerms) {
			totalEFunc.addTerm(term);
		}
		
		// low DoF single and pair (residues 1 and 13 and is ALA)
		MultiTermEnergyFunction alaEFunc = makeEFunc(singleTerms.get(0), "ALA");
		MultiTermEnergyFunction alaAlaEFunc = makeEFunc(pairTerms.get(getPairIndex(0, 12)), "ALA", "ALA");
		
		// high DoF single and pair (residues 9 and 22 are ARG)
		MultiTermEnergyFunction argEFunc = makeEFunc(singleTerms.get(8), "ARG");
		MultiTermEnergyFunction argArgEFunc = makeEFunc(pairTerms.get(getPairIndex(8, 21)), "ARG", "ARG");
		
		// BENCHMARKING: each test is run for a certain number of iterations
		// on Jeff's dual-core laptop (with hyper-threading, 2 physical cores, 4 logical cores),
		// the initial iterations were set to target about 20 seconds of running time
		// Jeff's laptop used the OpenJDK v1.8.0_72 JVM for these tests
		// run times should be at least 10 seconds so the JVM can "settle into" its internal optimizations
		// also, multi-threading was only used on the multi-term energy functions
		// it doesn't make much sense on the single term tests =P
		// we should expect a parallel speedup of close to 2x on this machine for 2 threads
		// for 4 threads, we shouldn't expect much higher than 2x speedup because the extra two "cores"
		// still share physical resources (like cache). I suspect we're mostly memory-bound here rather
		// than cpu-bound, so I'd bet a machine with 4 physical cores would do much better on 4 threads.
		// having more cache should speed things up a lot! =)
		
		// notation below (trialN values are operations per second):
		// nameOfEFunc: numThreads x numIters = [trial1, trial2, trial2]
		
		// BEFORE OPTIMIZATIONS (2016-04-13):
		// total: 1 x 40   = [2.05, 2.05, 2.04]
		// total: 4 x 40   = [3.63, 3.28, 3.51] (roughly 1.7x speedup... not great for 4 threads on 2 cores)
		// ala:   1 x 10e6 = [695400.29, 681174.63, 676398.07]
		// ala^2: 1 x 5e6  = [283561.50, 273957.20, 277309.58]
		// arg:   1 x 2e6  = [90464.21, 90977.86, 82896.86]
		// arg^2: 1 x 13e5 = [62300.29, 63856.59, 63150.76]
		
		// ugh... looks like results are inconsistent from day to day
		// that means I have to re-run the benchmarks every day  ;_;
		
		// 2016-04-14: benchmark, no optimizations (why is it so much slower than yesterday?!):
		// total: 1 x 40   = [1.70, 1.76, 1.68] avg=1.69
		// ala:   1 x 10e6 = [579805.52, 587120.14, 581297.23] avg=582740.95
		// ala^2: 1 x 5e6  = [236205.30, 232150.31, 233085.37] avg=233813.64
		// arg:   1 x 2e6  = [74910.83, 75621.81, 76753.69] avg=75762.08
		// arg^2: 1 x 13e5 = [54733.84, 54349.87, 53882.85] avg=54322.19
		
		// 2016-04-14: some work on the simplest test case...
		// no optimizations:
		// ala:   1 x 10e6 = [554713.03, 547208.29, 584274.57]
		// get rid of array allocations:
		// ala:   1 x 10e6 = [597278.34, 584392.16, 586566.90] => maybe slight improvement? hard to tell
		// inlining energy functions
		// ala:   1 x 10e6 = [591555.65, 570630.89, 591634.03] => no significant change
		// better handling of hydrogen electrostatics/vdW flags
		// ala:   1 x 10e6 = [645581.18, 618122.16, 621685.62] => modest but significant speedup
		// flatten solvation calculations to 1d and premultiply as much as possible
		// ala:   1 x 10e6 = [681132.80, 698386.85, 684975.58] => pretty noticeable speedup, not bad
		
		// 2016-04-14: after all the single-threaded optimizations I can think to do:
		// total: 1 x 40   = [1.77, 1.72, 1.74] avg=1.74, speedup=1.03x
		// ala:   1 x 10e6 = [707171.08, 680953.47, 699251.56] avg=695792.04, speedup=1.19x
		// ala^2: 1 x 5e6  = [253889.23, 255342.68, 268118.92] avg=259116.94, speedup=1.11x
		// arg:   1 x 2e6  = [84477.03, 81017.39, 82983.05] avg=82825.82, speedup=1.09x
		// arg^2: 1 x 13e5 = [56371.44, 53888.25, 54721.48] avg=54993.72, speedup=1.01x
		
		// 2016-04-15: working on parallelism, today's benchmarks pre-optimization:
		// total: 1 x 40   = [1.82, 1.80, 1.82]
		// total: 2 x 40   = [3.31, 3.34, 3.33] => 1.84x speedup, not bad! =)
		// total: 4 x 40   = [3.52, 3.55, 3.54] => 1.95x speedup, not great =(
		
		// 2016-04-15: custom parallel processors, before any thread-local memory optimizations:
		// total: 2 x 40   = [3.53, 3.56, 3.53] => 1.96x speedup, yeah! =D
		// total: 4 x 40   = [3.72, 3.73, 3.70] => 2.05x speedup, beats java's generic parallelism by a bit, makes me happy =)
		
		// 2016-04-15: ugh... something changed again. Need to redo the benchmarks. What is going on the with the JVM...
		// total: 1 x 40   = [1.88, 1.86, 1.87]
		// total: 2 x 40   = [3.58, 3.58, 3.58] => 1.91x speedup
		// total: 4 x 40   = [3.77, 3.78, 3.80] => 2.02x speedup
		
		// 2016-04-15: final post-optimization
		// total: 2 x 40   = [3.66, 3.63, 3.65] => 1.95x speedup
		// total: 4 x 40   = [3.81, 3.83, 3.84] => 2.05x speedup
		
		// 2016-04-18: today's benchmarks, all current optimizations:
		// total: 1 x 40  = [1.87, 1.87, 1.86]
		// total: 4 x 40  = [3.82, 3.78, 3.83]
		
		// 2016-04-18: added energy caching... performance is through the roof! =D
		// the caches are so fast, we're spending more time in thread synchronization now
		// total: 1 x 400  = [43.91, 41.82, 43.36] => 23.05x speedup over benchmark
		// total: 2 x 400  = [72.19, 75.89, 78.59] => 1.76x over single
		// total: 4 x 400  = [96.51, 95.59, 94.88] => 2.22x over single)
		
		// DO EEEEEETTT!!!
		final int thou = 1000;
		profile(totalEFunc, 400, -2937.319131349481300);
		//profile(alaEFunc, 10*thou*thou, -11.682132279211443);
		//profile(alaAlaEFunc, 5*thou*thou, 0.005174712669362);
		//profile(argEFunc, 2*thou*thou, -24.665020813395530);
		//profile(argArgEFunc, 1300*thou, 0.082934569827485);
	}
	
	private static int getPairIndex(int i, int j) {
		assert (i != j);
		if (i > j) {
			return i*(i - 1)/2 + j;
		} else {
			return j*(j - 1)/2 + i;
		}
	}
	
	private static MultiTermEnergyFunction makeEFunc(SingleResEnergy term, String expectedName) {
		checkNames(expectedName, term.getRes().template.name);
		System.out.println(String.format("%3s:     %5d terms",
			expectedName, term.getFFEnergy().getNumTerms()
		));
		return makeEFunc(term);
	}
	
	private static MultiTermEnergyFunction makeEFunc(ResPairEnergy term, String expectedName1, String expectedName2) {
		checkNames(
			expectedName1, term.getRes1().template.name,
			expectedName2, term.getRes2().template.name
		);
		System.out.println(String.format("%3s-%3s: %5d terms",
			expectedName1, expectedName2, term.getFFEnergy().getNumTerms()
		));
		return makeEFunc(term);
	}
	
	private static MultiTermEnergyFunction makeEFunc(EnergyFunction term) {
		MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
		efunc.addTerm(term);
		return efunc;
	}
	
	private static void checkNames(String ... names) {
		boolean isGood = true;
		for (int i=0; i<names.length/2; i++) {
			String expected = names[i*2 + 0];
			String observed = names[i*2 + 1];
			if (!expected.toLowerCase().equals(observed.toLowerCase())) {
				isGood = false;
				break;
			}
		}
		if (!isGood) {
			StringBuilder expected = new StringBuilder();
			StringBuilder observed = new StringBuilder();
			for (int i=0; i<names.length/2; i++) {
				if (i > 0) {
					expected.append(", ");
					observed.append(", ");
				}
				expected.append(names[i*2 + 0].toUpperCase());
				observed.append(names[i*2 + 1].toUpperCase());
			}
			throw new Error(String.format("Residue names don't match\n\texpected: %s\n\tobserved: %s", expected, observed));
		}
	}
	
	private static void profile(EnergyFunction efunc, int numIterations, double expectedEnergy) {

		System.out.println("\n\nStarting energy calculations...");
		
		// make sure the energy is correct before running the test
		checkEnergy(expectedEnergy, efunc.getEnergy());

		// time the energy calculations
		double energy = 0;
		Stopwatch stopwatch = new Stopwatch();
		stopwatch.start();
		for (int i=0; i<numIterations; i++) {
			energy = efunc.getEnergy();
			
			// DEBUG
			//checkEnergy(expectedEnergy, energy);
		}
		stopwatch.stop();
		
		System.out.println(String.format("Calculated %d energies in %s\n\te: %.15f\n\tOpS: %.2f",
			numIterations, stopwatch.getTime(TimeUnit.MILLISECONDS), energy,
			numIterations/stopwatch.getTimeS()
		));
		
		// make sure the energy is correct after the test too
		checkEnergy(expectedEnergy, energy);
	}
	
	private static void checkEnergy(double expected, double observed) {
		final double Epsilon = 1e-14;
		double relErr = Math.abs(expected - observed)/expected;
		if (relErr > Epsilon) {
			throw new Error(String.format(
				"Energy is wrong, undo that 'optimization'!\n\texpected: %.15f\n\tcalculated: %.15f\n\trelErr: %.15f",
				expected, observed, relErr
			));
		} else {
			System.out.println(String.format("Energy is correct\n\trelErr: %.15f", relErr));
		}
	}
}
