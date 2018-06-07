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

package edu.duke.cs.osprey.control;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSAllowedSeqs;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.impl.KSImplKAStar;
import edu.duke.cs.osprey.kstar.impl.KSImplLinear;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFFactory;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;


/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 * Deprecated, use the new KStar class instead
 *
 * Also, we're not maintaining or testing this code anymore,
 * and there's some evidence from failed tests (tests that have since been
 * removed) that it may no longer work correctly. Continue using this code
 * only with extreme skepticism.
 */
@Deprecated
public class KStarCalculator {

	KSConfigFileParser cfp;

	double Ew;//energy window for enumerating conformations: 0 for just GMEC
	double I0 = 0;//initial value of iMinDEE pruning interval
	boolean doIMinDEE;
	boolean useContFlex;

	HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs = new HashMap<>();
	
	
	public KStarCalculator( KSConfigFileParser cfgP ) {
		cfp = cfgP;

		Ew = cfp.params.getDouble("Ew", 5);
		doIMinDEE = cfp.params.getBool("imindee",false);
		if(doIMinDEE){
			I0 = cfp.params.getDouble("Ival", 5);
		}
		useContFlex = cfp.params.getBool("doMinimize",false);
		if(doIMinDEE && !useContFlex)
			throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility. "
					+ "Change the value of doMinimize to 'true'.");

		PositionConfSpace.dihedFlexInterval = cfp.params.getDouble("dihedFlexInterval");

		PFAbstract.suppressOutput = cfp.params.getBool("kStarPFuncSuppressOutput");
		PFAbstract.targetEpsilon = cfp.params.getDouble("epsilon");
		PFAbstract.setPhase2Method(cfp.params.getValue("kStarPhase2Method"));
		PFAbstract.qCapacity = cfp.params.getInt("kStarPFuncQCap", ThreadParallelism.getNumThreads()*8);

		PFAbstract.setCFGImpl(cfp.params.getValue("kStarPFuncMethod"));
		PFAbstract.setStabilityThresh( cfp.params.getDouble("kStarPFuncStabThresh") );
		PFAbstract.setConfsThreadBuffer( cfp.params.getInt("kStarPFuncConfsThreadBuffer", 4) );
		PFAbstract.setNumThreads( cfp.params.getInt("kStarPFuncThreads") );

		PFAbstract.saveTopConfsAsPDB = cfp.params.getBool("kStarSaveTopConfsAsPDB");
		PFAbstract.setNumTopConfsToSave( cfp.params.getInt("kStarNumTopConfsToSave") );
		PFAbstract.useMaxKSConfs = cfp.params.getBool("kStarUseMaxKSConfs");
		PFAbstract.setMaxKSconfs( cfp.params.getInt("kStarMaxKSConfs") );

		PFAbstract.setHotMethod( "kStarPFunctHotMethod", cfp.params.getValue("kStarPFunctHotMethod") );

		// check hots for validity
		if(!PFAbstract.getHotMethod().equalsIgnoreCase("none"))
			cfp.getHighOrderTuplesByPDBResNum();

		PFAbstract.setHotNumRes( "kStarPFuncHotNumRes", cfp.params.getInt("kStarPFuncHotNumRes", 3) );
		PFAbstract.setHotBoundPct( "kStarPFuncHotBoundPct", cfp.params.getDouble("kStarPFuncHotBoundPct", 0.03) );
		PFAbstract.setHotTopRotsPct( "KStarPFuncHotTopRotsPct", cfp.params.getDouble("kStarPFuncHotTopRotsPct", 0.0) );

		KSAbstract.runTimeout = cfp.params.getInt("kStarRunTimeout");
		KSAbstract.doCheckPoint = cfp.params.getBool("kStarDoCheckpoint");
		KSAbstract.setCheckPointInterval(cfp.params.getInt("kStarCheckpointInterval"));
		//KSAbstract.interMutationConst = cfp.params.getDouble("kStarInterMutationConst", 0.0);

		KSImplKAStar.useTightBounds = cfp.params.getBool("kStarUseTightBounds", true);
		KSImplKAStar.nodeExpansionMethod = cfp.params.getValue("kStarNodeExpansion", "parallel1");
	}


	protected ArrayList<String> getWTSequence() {
		return cfp.getWTSequence();
	}


	protected ArrayList<ArrayList<String>> getMutationsFromFile( String path ) throws IOException {

		if( !new File(path).exists() )
			throw new RuntimeException("ERROR: " + path + " does not exist");

		ArrayList<ArrayList<String>> ans = new ArrayList<>();
		try (BufferedReader br = new BufferedReader(new FileReader(path))) {
			String line;
			while ((line = br.readLine()) != null) {
				ArrayList<String> l = new ArrayList<String>();

				/*
				int pos = StringParsing.ordinalIndexOf(line, " ", 1) + 1;
				for( String s : Arrays.asList( line.substring(pos).split(" ") ) ) {
					l.add(s.trim());
				}
				 */

				line = line.trim();
				for( String s : Arrays.asList( line.split(" ") ) ) {
					if(s.contains("-")) s = s.split("-")[0];
					l.add(s.trim());
				}

				ans.add(l);
			}
		}

		if(ans.size() > 0) {
			// check for correct length
			KSAllowedSeqs pl = strand2AllowedSeqs.get(2);

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


	public ArrayList<ArrayList<String>> truncateAllowedSequences(String path) throws IOException {

		// read .mut file
		// filter list of mutations; only run those listed
		ArrayList<ArrayList<String>> mutations = getMutationsFromFile( path );

		if(mutations == null) 
			return null;

		KSAllowedSeqs pl = strand2AllowedSeqs.get(2);
		KSAllowedSeqs p = strand2AllowedSeqs.get(0);
		KSAllowedSeqs l = strand2AllowedSeqs.get(1);

		int plLen = pl.getSequenceLength(), pLen = p.getSequenceLength();

		pl.getStrandSeqList().clear(); 
		if(pl.addWT) pl.getStrandSeqList().add(pl.getWTSeq());
		
		p.getStrandSeqList().clear(); 
		if(p.addWT) p.getStrandSeqList().add(p.getWTSeq());
		
		l.getStrandSeqList().clear(); 
		if(l.addWT) l.getStrandSeqList().add(l.getWTSeq());

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

		KSAllowedSeqs complexSeqs = cfp.getAllowedSequences(2, null);
		strand2AllowedSeqs.put(2, complexSeqs);
		strand2AllowedSeqs.put(0, cfp.getAllowedSequences(0, complexSeqs));
		strand2AllowedSeqs.put(1, cfp.getAllowedSequences(1, complexSeqs));
	}

	private KSAbstract makeKStar() {
		switch (cfp.params.getValue("kStarMethod")) {
			case "kastar": return new KSImplKAStar(cfp);
			case "linear": return new KSImplLinear(cfp);
			default: throw new UnsupportedOperationException("ERROR: currently supported implementations are 'linear' and 'kastar'");
		}
	}

	public KSAbstract calcKStarScores() {
		
		// use the strand cache to reduce minimization overhead
		PFFactory.initStrandInfoCache();

		try {
			
			cfp.verifyStrandsMutuallyExclusive();

			generateAllowedSequences();

			String mutFilePath = cfp.params.getValue("mutfile", "");
			if(mutFilePath.length() > 0) {
				truncateAllowedSequences(mutFilePath);
			}
		} catch (IOException ex) {
			
			// don't make the caller trap this exception if it doesn't care,
			// but still give the caller a choice by wrapping the Exception with Error instead of calling System.exit()
			throw new Error(ex);
		}

		// run K*
		KSAbstract kstar = makeKStar();
		kstar.init(strand2AllowedSeqs);
		kstar.run();
		
		// cleanup strand cache after all pfuncs computed
		PFFactory.cleanupStrandInfoCache();
		
		return kstar;
	}

}
