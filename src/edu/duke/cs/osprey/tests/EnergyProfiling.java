package edu.duke.cs.osprey.tests;

import java.io.File;
import java.util.ArrayList;

import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;

public class EnergyProfiling {
	
	private static final long NSpS = 1000000000;
	private static final long NSpMS = 1000000;

	@SuppressWarnings("unused")
	public static void main(String[] args) {
		
		// check the cwd
		String path = new File("").getAbsolutePath();
		if (!path.endsWith("test/DAGK")) {
			throw new Error("This profiler was designed to run in the test/DAGK folder\n\tcwd: " + path);
		}

		// load configuration
		ConfigFileParser cfp = new ConfigFileParser(new String[] {"-c", "KStar.cfg"});
		cfp.loadData();
		
		// configure energy function parallelization
		MultiTermEnergyFunction.setNumThreads(1);
		if (MultiTermEnergyFunction.getNumThreads() > 1) {
			System.setProperty(
				"java.util.concurrent.ForkJoinPool.common.parallelism",
				Integer.toString(MultiTermEnergyFunction.getNumThreads())
			);
		}
		
		// read a big test protein, the bigger the better
		Molecule m = PDBFileReader.readPDBFile("2KDC.P.forOsprey.pdb");
		
		System.out.println("\n\nBuilding energy functions...");
		
		// get all the energy function terms
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		ArrayList<SingleResEnergy> singleTerms = new ArrayList<>();
		for (Residue res : m.residues) {
			singleTerms.add(new SingleResEnergy(res, ffparams));
		}
		ArrayList<ResPairEnergy> pairTerms = new ArrayList<>();
		for (int i=0; i<m.residues.size(); i++) {
			for (int j=0; j<i; j++) {
				pairTerms.add(new ResPairEnergy(m.residues.get(i), m.residues.get(j), ffparams));
			}
		}
		
		System.out.println(String.format("Built %d energy terms:\n\t%d singles\n\t%d pairs",
			singleTerms.size() + pairTerms.size(),
			singleTerms.size(),
			pairTerms.size()
		));

		// build a few different energy functions to test
		
		// the shole shebang
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
		// on Jeff's quad-core laptop, the initial iterations were set to target about 20 seconds of running time
		// run times should be at least 10 seconds so the JVM can "settle into" its internal optimizations
		// also, multi-threading was only used on the multi-term energy functions
		// it doesn't make much sense on the single term tests =P
		
		// notation below (trialN values are operations per second):
		// nameOfEFunc: numThreads x numIters = [trial1, trial2, trial2]
		
		// BEFORE OPTIMIZATIONS (2016-04-13):
		// total: 1 x 40   = [2.05, 2.05, 2.04]
		// total: 4 x 40   = [3.63, 3.28, 3.51] (roughly 1.6-1.7x speedup... not great for 4 threads=cores)
		// ala:   1 x 10e6 = [695400.29, 681174.63, 676398.07]
		// ala^2: 1 x 5e6  = [283561.50, 273957.20, 277309.58]
		// arg:   1 x 2e6  = [90464.21, 90977.86, 82896.86]
		// arg^2: 1 x 13e5 = [62300.29, 63856.59, 63150.76]
		
		// ugh... looks like results are inconsistent from day to day
		// that means I have to re-run the benchmarks every day  ;_;
		
		// 2016-04-13: benchmark, no optimizations:
		// total: 1 x 40   = [1.70, 1.76, 1.68] avg=1.69
		// ala:   1 x 10e6 = [579805.52, 587120.14, 581297.23] avg=582740.95
		// ala^2: 1 x 5e6  = [236205.30, 232150.31, 233085.37] avg=233813.64
		// arg:   1 x 2e6  = [74910.83, 75621.81, 76753.69] avg=75762.08
		// arg^2: 1 x 13e5 = [54733.84, 54349.87, 53882.85] avg=54322.19
		
		// 2016-04-13: some work on the simplest test case...
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
		
		// 2016-04-13: after all the single-threaded optimizations I can think to do:
		// total: 1 x 40   = [1.77, 1.72, 1.74] avg=1.74, speedup=1.03x
		// ala:   1 x 10e6 = [707171.08, 680953.47, 699251.56] avg=695792.04, speedup=1.19x
		// ala^2: 1 x 5e6  = [253889.23, 255342.68, 268118.92] avg=259116.94, speedup=1.11x
		// arg:   1 x 2e6  = [84477.03, 81017.39, 82983.05] avg=82825.82, speedup=1.09x
		// arg^2: 1 x 13e5 = [56371.44, 53888.25, 54721.48] avg=54993.72, speedup=1.01x
		
		// DO EEEEEETTT!!!
		final int thou = 1000;
		//profile(totalEFunc, 40, -2937.319131349481300);
		//profile(alaEFunc, 10*thou*thou, -11.682132279211443);
		//profile(alaAlaEFunc, 5*thou*thou, 0.005174712669362);
		//profile(argEFunc, 2*thou*thou, -24.665020813395530);
		profile(argArgEFunc, 1300*thou, 0.082934569827485);
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
		
		// make sure the energy is correct
		final double Epsilon = 1e-14;
		double energy = efunc.getEnergy();
		double err = Math.abs(expectedEnergy - energy);
		if (err > Epsilon) {
			throw new Error(String.format(
				"Energy is wrong, undo that 'optimization'!\n\texpected: %.15f\n\tcalculated: %.15f\n\terr: %.15f",
				expectedEnergy, energy, err
			));
		} else {
			System.out.println(String.format("Energy is correct\n\terr: %.15f", err));
		}

		// time the energy calculations
		long startNs = System.nanoTime();
		for (int i=0; i<numIterations; i++) {
			energy = efunc.getEnergy();
		}
		long diffNs = System.nanoTime() - startNs;
		
		System.out.println(String.format("Calculated %d energies in %d ms\n\te: %.15f\n\tOpS: %.2f",
			numIterations, diffNs/NSpMS, energy,
			(double)numIterations/diffNs*NSpS
		));
	}
}
