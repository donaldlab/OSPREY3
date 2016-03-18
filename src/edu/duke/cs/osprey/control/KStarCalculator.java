package edu.duke.cs.osprey.control;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import edu.duke.cs.osprey.kstar.AllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.Strand;
import edu.duke.cs.osprey.kstar.impl.KSImplLinear;
import edu.duke.cs.osprey.kstar.impl.KSImplAStar;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.minimization.MinimizerFactory;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
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

	HashMap<Integer, AllowedSeqs> strand2AllowedSeqs = new HashMap<>();

	public KStarCalculator ( ConfigFileParser cfgP ) {
		cfp = cfgP;

		Ew = cfp.getParams().getDouble("Ew", 5);
		doIMinDEE = cfp.getParams().getBool("imindee",false);
		if(doIMinDEE){
			I0 = cfp.getParams().getDouble("Ival", 5);
		}
		useContFlex = cfp.getParams().getBool("doMinimize",false);
		if(doIMinDEE && !useContFlex)
			throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility. "
					+ "Change the value of doMinimize to 'true'.");
		
		PFAbstract.targetEpsilon = cfp.getParams().getDouble("epsilon", 0.03);
		PFAbstract.qCapacity = cfp.getParams().getInt("pFuncQCap", 4194304);
		PFAbstract.waitUntilCapacity = cfp.getParams().getBool("pFuncQWait", false);

		PFAbstract.eMinMethod = cfp.getParams().getValue("eMinMethod", "ccd");
		PFAbstract.setImpl(cfp.getParams().getValue("pFuncMethod", "1npcpmcache"));
		PFAbstract.setStabilityThresh( cfp.getParams().getDouble("pFuncStabThresh", 0.00001) );
		PFAbstract.setConfsThreadBuffer( cfp.getParams().getInt("pFuncConfsThreadBuffer", 4) );
		PFAbstract.setNumFibers( cfp.getParams().getInt("pFuncFibers", 1) );
		PFAbstract.setNumThreads( cfp.getParams().getInt("pFuncThreads", ThreadParallelism.getNumThreads()) );
		PFAbstract.setServerList( cfp.getParams().getValue("pFuncServerList", "localhost").split("\\s+") );
		PFAbstract.setNumRemoteClients( cfp.getParams().getInt("pFuncClients", 1) );

		PFAbstract.saveTopConfsAsPDB = cfp.getParams().getBool("saveTopConfsAsPDB", false);
		PFAbstract.setNumTopConfsToSave( cfp.getParams().getInt("numTopConfsToSave", 10) );
		PFAbstract.useMaxKSConfs = cfp.getParams().getBool( "useMaxKSConfs", false );
		PFAbstract.setMaxKSconfs( cfp.getParams().getInt("maxKSconfs", 100000) );

		MinimizerFactory.setImpl( PFAbstract.eMinMethod );
		
		KSAbstract.refinePInterval = cfp.getParams().getBool("KStarRefinePruning", false);
		KSAbstract.doCheckpoint = cfp.getParams().getBool("doKStarCheckpoint", false);
		KSAbstract.setCheckPointInterval(cfp.getParams().getInt("KStarCheckpoint", 100000));
	}


	protected ArrayList<String> getWTSequence() {
		return cfp.getWTSequence();
	}


	protected ArrayList<ArrayList<String>> readMutFile( String path ) throws Exception {

		if( !new File(path).exists() )
			throw new RuntimeException("ERROR: " + path + " does not exist");

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


	public ArrayList<ArrayList<String>> filterByMutFile(String path) throws Exception {

		// read .mut file
		// filter list of mutations; only run those listed
		ArrayList<ArrayList<String>> mutations = readMutFile( path );

		if(mutations == null) 
			return null;

		AllowedSeqs pl = strand2AllowedSeqs.get(Strand.COMPLEX);
		AllowedSeqs p = strand2AllowedSeqs.get(Strand.PROTEIN);
		AllowedSeqs l = strand2AllowedSeqs.get(Strand.LIGAND);

		int plLen = pl.getSequenceLength(), pLen = p.getSequenceLength();

		if(mutations.get(0).size() != plLen) {
			throw new RuntimeException("ERROR: mutfile sequences have length " + mutations.get(0).size()
					+ " but sequences in this design have length " + plLen);
		}

		ArrayList<String> plWT = pl.getStrandSeqAtPos(0);
		ArrayList<String> pWT = p.getStrandSeqAtPos(0);
		ArrayList<String> lWT = l.getStrandSeqAtPos(0);

		pl.getStrandSeqList().clear(); pl.getStrandSeqList().add(plWT);
		p.getStrandSeqList().clear(); p.getStrandSeqList().add(pWT);
		l.getStrandSeqList().clear(); l.getStrandSeqList().add(lWT);

		for(ArrayList<String> list : mutations) {
			
			ArrayList<String> seq = AllowedSeqs.addPosToSeq(list, pl.getFlexRes());
			
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

		pl.truncateAllowedAAs();
		p.truncateAllowedAAs();
		l.truncateAllowedAAs();

		return mutations;
	}


	private void createAllowedSequences() {
		AllowedSeqs complexSeqs = cfp.getAllowedSequences(Strand.COMPLEX, null);
		strand2AllowedSeqs.put(Strand.COMPLEX, complexSeqs);
		strand2AllowedSeqs.put(Strand.PROTEIN, cfp.getAllowedSequences(Strand.PROTEIN, complexSeqs));
		strand2AllowedSeqs.put(Strand.LIGAND, cfp.getAllowedSequences(Strand.LIGAND, complexSeqs));
	}


	public void calcKStarScores() {

		try {

			createAllowedSequences();

			String mutFilePath = cfp.getParams().getValue("mutfile", "");
			if(mutFilePath.length() > 0) {
				filterByMutFile(mutFilePath);
			}

			String ksMethod = cfp.getParams().getValue("kstarmethod", "linear");

			switch( ksMethod ) {

			case "astar":
				KSImplAStar astar = new KSImplAStar(cfp);
				astar.init(strand2AllowedSeqs);
				astar.run();
				break;
			
			case "linear":
			default:
				KSImplLinear linear = new KSImplLinear(cfp);
				linear.init(strand2AllowedSeqs);
				linear.run();
				break;
			}

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

	}

}
