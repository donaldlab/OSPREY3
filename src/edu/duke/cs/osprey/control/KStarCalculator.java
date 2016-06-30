package edu.duke.cs.osprey.control;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import edu.duke.cs.osprey.kstar.KSAllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.Termini;
import edu.duke.cs.osprey.kstar.impl.KSImplLinear;
import edu.duke.cs.osprey.kstar.impl.KSImplKAStar;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFNew02;
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

	HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs = new HashMap<>();

	public KStarCalculator( ConfigFileParser cfgP ) {
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
		
		PFAbstract.suppressOutput = cfp.getParams().getBool("pFuncSuppressOutput", false);
		PFAbstract.targetEpsilon = cfp.getParams().getDouble("epsilon", 0.03);
		PFAbstract.setPhase2Method(cfp.getParams().getValue("kStarPhase2Method", "fast"));
		PFAbstract.qCapacity = cfp.getParams().getInt("pFuncQCap", ThreadParallelism.getNumThreads()*8);

		PFAbstract.eMinMethod = cfp.getParams().getValue("eMinMethod", "ccd");
		PFAbstract.setCFGImpl(cfp.getParams().getValue("pFuncMethod", new PFNew02().getImpl()));
		PFAbstract.setStabilityThresh( cfp.getParams().getDouble("pFuncStabThresh", 0) );
		PFAbstract.setConfsThreadBuffer( cfp.getParams().getInt("pFuncConfsThreadBuffer", 4) );
		PFAbstract.setNumFibers( cfp.getParams().getInt("pFuncFibers", 1) );
		PFAbstract.setNumThreads( cfp.getParams().getInt("pFuncThreads", ThreadParallelism.getNumThreads()) );
		PFAbstract.setServerList( cfp.getParams().getValue("pFuncServerList", "localhost").split("\\s+") );
		PFAbstract.setNumRemoteClients( cfp.getParams().getInt("pFuncClients", 1) );

		PFAbstract.saveTopConfsAsPDB = cfp.getParams().getBool("saveTopConfsAsPDB", false);
		PFAbstract.setNumTopConfsToSave( cfp.getParams().getInt("numTopConfsToSave", 10) );
		PFAbstract.useMaxKSConfs = cfp.getParams().getBool("useMaxKSConfs", false);
		PFAbstract.setMaxKSconfs( cfp.getParams().getInt("maxKSconfs", 100000) );
		
		PFAbstract.setHotMethod( "pFunctHotMethod", cfp.getParams().getValue("pFunctHotMethod", "none") );
		
		// check hots for validity
		if(!PFAbstract.getHotMethod().equalsIgnoreCase("none"))
			cfp.getHighOrderTuplesByPDBResNum();
		
		PFAbstract.setHotNumRes( "pFuncHotNumRes", cfp.getParams().getInt("pFuncHotNumRes", 3) );
		PFAbstract.setHotBoundPct( "pFuncHotBoundPct", cfp.getParams().getDouble("pFuncHotBoundPct", 0.03) );
		PFAbstract.setHotTopRotsPct( "pFuncHotTopRotsPct", cfp.getParams().getDouble("pFuncHotTopRotsPct", 0.0) );
		
		MinimizerFactory.setImpl( PFAbstract.eMinMethod );

		KSAbstract.doCheckPoint = cfp.getParams().getBool("doKStarCheckpoint", false);
		KSAbstract.setCheckPointInterval(cfp.getParams().getInt("kStarCheckpoint", 50000));
		//KSAbstract.interMutationConst = cfp.getParams().getDouble("kStarInterMutationConst", 0.0);
		
		KSImplKAStar.useTightBounds = cfp.getParams().getBool("kStarUseTightBounds", true);
		KSImplKAStar.nodeExpansionMethod = cfp.getParams().getValue("kStarNodeExpansion", "parallel_1");
	}


	protected ArrayList<String> getWTSequence() {
		return cfp.getWTSequence();
	}


	protected ArrayList<ArrayList<String>> getMutationsFromFile( String path ) throws Exception {

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

		if(ans.size() > 0) {
			// check for correct length
			KSAllowedSeqs pl = strand2AllowedSeqs.get(Termini.COMPLEX);
			
			for(ArrayList<String> seq : ans) {
				if(seq.size() != pl.getSequenceLength())
					throw new RuntimeException("ERROR: sequence " + KSAbstract.list1D2String(seq, " ") 
							+ " has a the wrong length");
			}
			
			// add residue numbers
			for(int i = 0; i < ans.size(); ++i) {
				ArrayList<String> seq = KSAllowedSeqs.addPosToSeq(ans.get(i), pl.getFlexRes());
				if(!pl.isAllowed(seq)) {
					throw new RuntimeException("ERROR: sequence " + KSAbstract.list1D2String(ans.get(i), " ") + 
							" is not allowed in the design space.\n Change resAllowed.");
				}
				ans.set(i, seq);
			}
			return ans;
		}
		
		return null;
	}


	public ArrayList<ArrayList<String>> truncateAllowedSequences(String path) throws Exception {

		// read .mut file
		// filter list of mutations; only run those listed
		ArrayList<ArrayList<String>> mutations = getMutationsFromFile( path );

		if(mutations == null) 
			return null;

		KSAllowedSeqs pl = strand2AllowedSeqs.get(Termini.COMPLEX);
		KSAllowedSeqs p = strand2AllowedSeqs.get(Termini.PROTEIN);
		KSAllowedSeqs l = strand2AllowedSeqs.get(Termini.LIGAND);

		int plLen = pl.getSequenceLength(), pLen = p.getSequenceLength();

		pl.getStrandSeqList().clear(); pl.getStrandSeqList().add(pl.getWTSeq());
		p.getStrandSeqList().clear(); p.getStrandSeqList().add(p.getWTSeq());
		l.getStrandSeqList().clear(); l.getStrandSeqList().add(l.getWTSeq());

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

		pl.truncateAllowedAAs();
		p.truncateAllowedAAs();
		l.truncateAllowedAAs();

		return mutations;
	}


	private void generateAllowedSequences() {
		
		KSAllowedSeqs complexSeqs = cfp.getAllowedSequences(Termini.COMPLEX, null);
		strand2AllowedSeqs.put(Termini.COMPLEX, complexSeqs);
		strand2AllowedSeqs.put(Termini.PROTEIN, cfp.getAllowedSequences(Termini.PROTEIN, complexSeqs));
		strand2AllowedSeqs.put(Termini.LIGAND, cfp.getAllowedSequences(Termini.LIGAND, complexSeqs));
	}


	public void calcKStarScores() {

		try {

			cfp.verifyStrandsMutuallyExclusive();
			
			generateAllowedSequences();

			String mutFilePath = cfp.getParams().getValue("mutfile", "");
			if(mutFilePath.length() > 0) {
				truncateAllowedSequences(mutFilePath);
			}

			String ksMethod = cfp.getParams().getValue("kstarmethod", "linear");

			switch( ksMethod ) {

			case "kastar":
				KSImplKAStar kastar = new KSImplKAStar(cfp);
				kastar.init(strand2AllowedSeqs);
				kastar.run();
				break;
			
			case "linear":
				KSImplLinear linear = new KSImplLinear(cfp);
				linear.init(strand2AllowedSeqs);
				linear.run();
				break;
				
			default:
				throw new RuntimeException("ERROR: currently supported implementations are 'linear' and 'kastar'");
			}

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

	}

}
