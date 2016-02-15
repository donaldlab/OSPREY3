package edu.duke.cs.osprey.control;


import java.util.ArrayList;
import java.util.HashMap;

import edu.duke.cs.osprey.confspace.AllowedSeqs;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.kstar.KSCalcManagerLinear;
import edu.duke.cs.osprey.kstar.PFAbstract;
import edu.duke.cs.osprey.minimization.MinimizerFactory;
import edu.duke.cs.osprey.pruning.PruningControl;


/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class KStarCalculator {

	ConfigFileParser cfp;

	double Ew;//energy window for enumerating conformations: 0 for just GMEC
	double I0 = 0;//initial value of iMinDEE pruning interval
	boolean doIMinDEE;
	boolean useContFlex;

	HashMap<Integer, SearchProblem> strand2SearchProblem = new HashMap<>();
	HashMap<Integer, AllowedSeqs> strand2AllowedSeqs = new HashMap<>();
	HashMap<Integer, PruningControl> strand2Pruning = new HashMap<>();

	public KStarCalculator ( ConfigFileParser cfgP ) {
		cfp = cfgP;

		Ew = cfp.params.getDouble("Ew", 0);
		doIMinDEE = cfp.params.getBool("imindee",false);
		if(doIMinDEE){
			I0 = cfp.params.getDouble("Ival",5);
		}
		useContFlex = cfp.params.getBool("doMinimize",false);
		if(doIMinDEE && !useContFlex)
			throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility");

		PFAbstract.targetEpsilon = cfp.params.getDouble("epsilon", 0.03);
		PFAbstract.rho = PFAbstract.targetEpsilon / (1.0 - PFAbstract.targetEpsilon);
		PFAbstract.qCapacity = cfp.params.getInt("pFuncQCap", (int)Math.pow(2, 17));
		PFAbstract.waitUntilCapacity = cfp.params.getBool("pFuncQWait", false);

		PFAbstract.eMinMethod = cfp.params.getValue("eMinMethod", "hbfgsccd");
		PFAbstract.setImplementation(cfp.params.getValue("pFuncMethod", "1npcpmcache"));
		PFAbstract.setStabilityThreshold( cfp.params.getDouble("pFuncStabilityThreshold", 1.0) );
		PFAbstract.setInterval( cfp.params.getValue("pFuncInterval", "0.01") );
		PFAbstract.setConfsThreadBuffer( cfp.params.getInt("pFuncConfsThreadBuffer", 4) );
		PFAbstract.setNumFibers( cfp.params.getInt("pFuncFibers", 4) );
		PFAbstract.setNumThreads( cfp.params.getInt("pFuncThreads", 2) );
		PFAbstract.setServerList( cfp.params.getValue("pFuncServerList", "localhost").split("\\s+") );
		PFAbstract.setNumRemoteClients( cfp.params.getInt("pFuncClients", 1) );

		MinimizerFactory.setImplementation( PFAbstract.eMinMethod );
	}


	private void createSearchProblems() {

		// computing allowed sequences can be time intensive for large sequence
		// spaces, so re-use it for different strands once it has been computed
		// for the complex.
		AllowedSeqs complexSeqs = null;
		strand2AllowedSeqs.put(Strand.COMPLEX, cfp.getAllowedSequences(Strand.COMPLEX, complexSeqs));
		complexSeqs = strand2AllowedSeqs.get(Strand.COMPLEX);
		strand2SearchProblem.put(Strand.COMPLEX, cfp.getSearchProblem(Strand.COMPLEX, complexSeqs));

		ArrayList<Integer> strands = new ArrayList<>();
		for( int strand = cfp.getNumStrands()-1; strand >= 0; --strand ) strands.add(strand);

		for( int strand : strands ) {
			strand2AllowedSeqs.put(strand, cfp.getAllowedSequences(strand, complexSeqs));

			System.out.println("\nCreating search problem for strand " + Strand.getStrandString(strand));
			strand2SearchProblem.put(strand, cfp.getSearchProblem(strand, strand2AllowedSeqs.get(strand)));
		}
	}


	private void createEnergyMatrices( boolean parallel ) {

		ArrayList<Integer> strands = new ArrayList<>();
		for( int strand = cfp.getNumStrands()-1; strand >= 0; --strand ) strands.add(strand);
		strands.add(Strand.COMPLEX);
		
		if( !parallel ) {
			// load sequentially
			for( int strand : strands ) {
				System.out.println("\nCreating energy matrix for strand " + Strand.getStrandString(strand));
				strand2SearchProblem.get(strand).loadEnergyMatrix();
			}
		}

		else {
			// load in parallel
			ArrayList<SearchProblem> sps = new ArrayList<>(); for( SearchProblem sp : strand2SearchProblem.values() ) sps.add(sp);
			ArrayList<Integer> indexes = new ArrayList<>(); for(int i = 0; i < sps.size(); ++i) indexes.add(i);
			indexes.parallelStream().forEach( i -> {
				SearchProblem sp = sps.get(i);
				sp.loadEnergyMatrix();
			});
		}
	}


	protected ArrayList<String> getWTSequence() {
		return cfp.getWTSequence();
	}


	private void pruneEnergyMatrices() {

		ArrayList<Integer> strands = new ArrayList<>();
		for( int strand = cfp.getNumStrands()-1; strand >= 0; --strand ) strands.add(strand);
		strands.add(Strand.COMPLEX);

		for( int strand : strands ) {
			System.out.println("\nPruning rotamers in strand " + Strand.getStrandString(strand));
			strand2Pruning.put(strand, cfp.getPruningControl(strand2SearchProblem.get(strand), Ew+I0, false, false));
			strand2Pruning.get(strand).prune();
		}
	}


	public void calcKStarScores() {

		try {

			createSearchProblems();
			//boolean parallel = true; 
			//createEnergyMatrices(parallel);
			//pruneEnergyMatrices();

			KSCalcManagerLinear kcm = null;
			String method = cfp.getParams().getValue("kstarmethod", "linear");
			switch( method ) {

			case "linear":
			default:
				kcm = new KSCalcManagerLinear(strand2SearchProblem, strand2AllowedSeqs, strand2Pruning, cfp);
				break;
			}

			kcm.run();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

	}

}
