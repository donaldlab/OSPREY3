package edu.duke.cs.osprey;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class EnergyMatrixProfiling {

	public static void main(String[] args)
	throws Exception {
		
		// check the cwd
		String path = new File("").getAbsolutePath();
		if (!path.endsWith("examples/DAGK")) {
			throw new Error("This profiler was designed to run in the examples/DAGK folder\n\tcwd: " + path);
		}

		// load configuration
		ConfigFileParser cfp = new ConfigFileParser();
		cfp.loadData();
		
		// init a conf space with lots of flexible residues, but no mutations
		final int NumFlexible = 20;
		ArrayList<String> flexRes = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<NumFlexible; i++) {
			flexRes.add(Integer.toString(i + 1));
			allowedAAs.add(new ArrayList<String>());
		}
		boolean addWt = true;
		boolean doMinimize = false;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = true;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"energyMatrixProfiling",
			"2KDC.P.forOsprey.pdb", 
			flexRes, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null, 
                        false, new ArrayList<>()
		);
		
		// compute the energy matrix
		System.out.println("\nComputing energy matrix...");
		Stopwatch stopwatch = new Stopwatch();
		stopwatch.start();
		EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
		emCalc.calcPEM();
		EnergyMatrix mat = emCalc.getEMatrix();
		stopwatch.stop();
		System.out.println("finished in " + stopwatch.getTime(TimeUnit.SECONDS, 2));
		
		/* oi... turns out writing to an encrypted folder is slow no matter what the code does
		   profiling this is just silly now... how about I don't write to an encrypted folder instead?
		// get a temp file for the energy matrix
		File tempFile = File.createTempFile("energyMatrixProfiling", "mat");
		try {
		
			System.out.println("\nWriting energy matrix...");
			Stopwatch.start();
			ObjectIO.writeObject(mat, tempFile.getAbsolutePath());
			Stopwatch.stop();
			System.out.println("finished in " + Stopwatch.getTime(TimeUnit.SECONDS, 2));
			
			System.out.println("\nReading energy matrix...");
			Stopwatch.start();
			ObjectIO.readObject(tempFile.getAbsolutePath(), true);
			Stopwatch.stop();
			System.out.println("finished in " + Stopwatch.getTime(TimeUnit.SECONDS, 2));
			
		} finally {
			
			// cleanup
			tempFile.delete();
		}
		*/
		
		// notation below (trialN values are operations per second):
		// numResidues,numIters: [trial1, trial2, trial2]
		
		// 2016-04-27
		// BEFORE OPTIMIZATIONS
		// 20,1e6: [27445.92, 28491.04, 27500.22]
		
		// THEORETICAL MAX (getXXX methods just return null)
		// 20,1e6: [211083.19, 222771.69, 223615.84] => about a 7.88x speedup
		
		// 2016-05-04
		// BEFORE OPTIMIZATIONS (re-benchmarking for today)
		// 20,1e6: [26702.41, 28462.70, 27122.50]
		
		// flatten tuple matrix to 1d
		// 20,1e6: [44531.32, 44333.52, 44617.89] => about a 1.62x speedup over benchmark
		
		// make double-specialized subclass
		// honestly, I didn't expect this to make much difference, but apparently it does
		// 20,1e6: [56037.08, 56238.11, 56288.57] => about a 2.05x speedup over benchmark
		
		// do lots of lookups in every spot
		System.out.println("\nProfiling reads...");
		stopwatch = new Stopwatch();
		stopwatch.start();
		final int n = 1000 * 1000;
		int numPos = mat.getNumPos();
		for (int i=0; i<n; i++) {
			for (int res1=0; res1<numPos; res1++) {
				int m1 = mat.getNumConfAtPos(res1);
				for (int i1=0; i1<m1; i1++) {
					mat.getOneBody(res1, i1);
					for (int res2=0; res2<res1; res2++) {
						int m2 = mat.getNumConfAtPos(res2);
						for (int i2=0; i2<m2; i2++) {
							mat.getPairwise(res1, i1, res2, i2);
						}
					}
				}
			}
		}
		stopwatch.stop();
		System.out.println("finished in " + stopwatch.getTime(TimeUnit.SECONDS, 2));
		System.out.println(String.format("\t%.2f OpS", n/stopwatch.getTimeS()));
	}
}
