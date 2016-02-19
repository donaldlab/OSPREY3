package edu.duke.cs.osprey.control;


import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import edu.duke.cs.osprey.confspace.AllowedSeqs;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.kstar.KSImplementationLinear;
import edu.duke.cs.osprey.kstar.PFAbstract;
import edu.duke.cs.osprey.minimization.MinimizerFactory;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.StringParsing;


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
		
		PFAbstract.saveTopConfsAsPDB = cfp.params.getBool("saveTopConfsAsPDB", false);
		PFAbstract.setNumTopConfsToSave( cfp.params.getInt("numTopConfsToSave", 10));
		
		MinimizerFactory.setImplementation( PFAbstract.eMinMethod );
	}


	private void createSearchProblems() {

		// computing allowed sequences can be time intensive for large sequence
		// spaces, so re-use it for different strands once it has been computed
		// for the complex.
		AllowedSeqs complexSeqs = cfp.getAllowedSequences(Strand.COMPLEX, null);
		strand2AllowedSeqs.put(Strand.COMPLEX, complexSeqs);
		System.out.println("\nCreating search problem for strand " + Strand.COMPLEX);
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


	protected ArrayList<ArrayList<String>> readMutFile( String path ) throws Exception {	
		ArrayList<ArrayList<String>> ans = new ArrayList<>();
		try (BufferedReader br = new BufferedReader(new FileReader(path))) {
			String line;
			while ((line = br.readLine()) != null) {
				ArrayList<String> l = new ArrayList<String>();

				int pos = StringParsing.ordinalIndexOf(line, " ", 1) + 1;
				for( String s : Arrays.asList( line.substring(pos).split(" ") ) ) {
					l.add(s.trim());
				}

				ans.add(l);
			}
		}

		if(ans.size() > 0)
			return ans;
		return null;
	}


	public void calcKStarScores() {

		try {
			createSearchProblems();
			
			// read .mut file
			// filter list of mutations; only run those listed
			ArrayList<ArrayList<String>> mutations = readMutFile( cfp.getParams().getValue("mutfile") );
			if(mutations != null) {
				
				AllowedSeqs pl = strand2AllowedSeqs.get(Strand.COMPLEX);
				AllowedSeqs p = strand2AllowedSeqs.get(Strand.PROTEIN);
				AllowedSeqs l = strand2AllowedSeqs.get(Strand.LIGAND);
				
				int plLen = pl.getSequenceLength(), pLen = p.getSequenceLength();
				
				if(mutations.get(0).size() != plLen) {
					throw new RuntimeException("ERROR: mutfile sequences have length " + mutations.get(0).size()
							+ " but sequences in this design have length " + plLen);
				}
				
				ArrayList<String> plWT = pl.getStrandSeq(0);
				ArrayList<String> pWT = p.getStrandSeq(0);
				ArrayList<String> lWT = l.getStrandSeq(0);
				
				pl.getStrandSeqList().clear(); pl.getStrandSeqList().add(plWT);
				p.getStrandSeqList().clear(); p.getStrandSeqList().add(pWT);
				l.getStrandSeqList().clear(); l.getStrandSeqList().add(lWT);

				for(ArrayList<String> seq : mutations) {
					if(!pl.getStrandSeqList().contains(seq)) {
						pl.getStrandSeqList().add(seq);
						
						// p
						ArrayList<String> pSubList = new ArrayList<>();
						for(String s : seq.subList(0, pLen)) pSubList.add(s);
						p.getStrandSeqList().add(pSubList);
						
						// l
						ArrayList<String> lSubList = new ArrayList<>();
						for(String s : seq.subList(pLen, plLen)) lSubList.add(s);
						l.getStrandSeqList().add(lSubList);
					}
				}
			}
			
			//boolean parallel = true; 
			//createEnergyMatrices(parallel);
			//pruneEnergyMatrices();

			KSImplementationLinear kcm = null;
			String method = cfp.getParams().getValue("kstarmethod", "linear");
			switch( method ) {

			case "linear":
			default:
				kcm = new KSImplementationLinear(mutations, strand2SearchProblem, strand2AllowedSeqs, strand2Pruning, cfp);
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
