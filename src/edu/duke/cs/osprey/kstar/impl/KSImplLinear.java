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

package edu.duke.cs.osprey.kstar.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import edu.duke.cs.osprey.kstar.KSAllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSCalc;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class KSImplLinear extends KSAbstract {
	
	private boolean contSCFlex;
	private String pfImpl;

	public KSImplLinear( KSConfigFileParser cfp ) {
		super( cfp );
	}

	public void init( HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs ) {

		this.strand2AllowedSeqs = strand2AllowedSeqs;
                checkAPPP();

		printSequences();

		createOutputDir();
		
		createEmatDir();

		if(doCheckPoint) createCheckPointDir();

		contSCFlex = cfp.params.getBool("doMinimize", true);
		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(contSCFlex));
		createEmats(contSCFlexVals);
		pfImpl = cfp.params.getValue("kStarPFuncMethod");
	}


	@Override
	public String getKSMethod() {
		return "linear";
	}


	@Override
	public void run() {

		if(strand2AllowedSeqs == null)
			throw new RuntimeException("ERROR: call init() method on this object before invoking run()");

		setStartTime(System.currentTimeMillis());

		if(doCheckPoint)
			runRR();

		else
			runFCFS();

		// print statistics
		System.out.println("\nK* calculations computed: " + getNumSeqsCompleted(0));
		System.out.println("Total # sequences: " + strand2AllowedSeqs.get(2).getNumSeqs());
		System.out.println("K* conformations processed: " + countProcessedConfs());
		System.out.println("Total # of conformations in search space: " + countTotNumConfs());
		System.out.println("K* running time: " + (System.currentTimeMillis()-getStartTime())/1000 + " seconds\n");
	}


	protected void runFCFS() {

		// each value corresponds to the desired flexibility of the 
		// pl, p, and l conformation spaces, respectively
		ArrayList<ArrayList<String>> strandSeqs = null;	
		ArrayList<Boolean> contSCFlexVals = new ArrayList<>(Arrays.asList(contSCFlex, contSCFlex, contSCFlex));
		ArrayList<String> pfImplVals = new ArrayList<>(Arrays.asList(pfImpl, pfImpl, pfImpl));

		// run wt
		int startSeq = 0;
		if( doWTCalc ) {
			computeWTCalc();
			startSeq = 1;
		}
			
		int numSeqs = strand2AllowedSeqs.get(2).getNumSeqs();
		for( int i = startSeq; i < numSeqs; ++i ) {

			// wt is seq 0, mutants are others
			System.out.println("\nComputing K* for sequence " + i + "/" + 
					(numSeqs-1) + ": " + 
					list1D2String(strand2AllowedSeqs.get(2).getStrandSeqAtPos(i), " ") + "\n");

			// get sequences
			strandSeqs = getStrandStringsAtPos(i);

			// create partition functions
			ConcurrentHashMap<Integer, PFAbstract> pfs = createPFs4Seqs(strandSeqs, contSCFlexVals, pfImplVals);

			// create K* calculation for sequence
			KSCalc calc = new KSCalc(i, pfs);

			// compute partition functions
			calc.run(wtKSCalc, false, true);

			// compute K* scores and print output if all 
			// partition functions are computed to epsilon accuracy
			if( calc.getEpsilonStatus() == EApproxReached.TRUE || calc.getEpsilonStatus() == EApproxReached.NOT_POSSIBLE ) {
				calc.printSummary( getOputputFilePath(), getStartTime(), getNumSeqsCreated(1), getNumSeqsCompleted(1) );
			}
		}
	}


	protected void runRR() {

		// each value corresponds to the desired flexibility of the 
		// pl, p, and l conformation spaces, respectively
		ArrayList<ArrayList<String>> strandSeqs = null;	
		ArrayList<Boolean> contSCFlexVals = new ArrayList<>(Arrays.asList(contSCFlex, contSCFlex, contSCFlex));
		ArrayList<String> pfImplVals = new ArrayList<>(Arrays.asList(pfImpl, pfImpl, pfImpl));

		// get all sequences		
		@SuppressWarnings("unchecked")
		HashSet<ArrayList<String>> seqSet = new HashSet<>((ArrayList<ArrayList<String>>) 
				ObjectIO.deepCopy(strand2AllowedSeqs.get(2).getStrandSeqList()));

		// run wt
		if( doWTCalc ) computeWTCalc();
		
		// remove completed seqs from the set of calculations we must compute
		seqSet.removeAll(getSeqsFromFile(getOputputFilePath()));

		int numSeqs = strand2AllowedSeqs.get(2).getNumSeqs();

		do {

			for( int i = 0; i < numSeqs; ++i ) {

				// wt is seq 0, mutants are others
				ArrayList<String> seq = strand2AllowedSeqs.get(2).getStrandSeqAtPos(i);

				if( !seqSet.contains(seq) ) continue;

				System.out.println("\nResuming K* for sequence " + i + "/" + (numSeqs-1) + ": " + list1D2String(seq, " ") + "\n");

				// get sequences
				strandSeqs = getStrandStringsAtPos(i);

				// create partition functions
				ConcurrentHashMap<Integer, PFAbstract> pfs = createPFs4Seqs(strandSeqs, contSCFlexVals, pfImplVals);

				// create K* calculation for sequence
				KSCalc calc = new KSCalc(i, pfs);

				// compute partition functions
				calc.run(wtKSCalc, false, true);

				// serialize P and L; should only happen once
				for( int strand : Arrays.asList(1, 0) ) {
					PFAbstract pf = calc.getPF(strand);
					if(!pf.checkPointExists()) calc.serializePF(strand);
				}

				// remove entry from checkpoint file
				PFAbstract pf = calc.getPF(2);
				calc.deleteSeqFromFile( pf.getSequence(), getCheckPointFilePath() );

				if( calc.getEpsilonStatus() != EApproxReached.FALSE ) {
					// we are not going to checkpoint this, so clean up
					calc.deleteCheckPointFile(2);
					seqSet.remove(seq);

					if( calc.getEpsilonStatus() == EApproxReached.TRUE || calc.getEpsilonStatus() == EApproxReached.NOT_POSSIBLE ) {
						calc.printSummary( getOputputFilePath(), getStartTime(), getNumSeqsCreated(1), getNumSeqsCompleted(1) );
					}
				}

				else {
					// remove partition funtion from memory, write checkpoint
					name2PF.remove(pf.getReducedSearchProblemName());
					calc.serializePF(2);
					calc.printSummary( getCheckPointFilePath(), getStartTime(), getNumSeqsCreated(0), getNumSeqsCompleted(0) );
				}
			}

		} while(seqSet.size() > 0);

		// delete checkpoint dir and checkpoint file
		ObjectIO.delete(getCheckPointDir());
	}

	public PFAbstract getPartitionFunction(int strand, ArrayList<String> seq) {
		return name2PF.get(getSearchProblemName(contSCFlex, strand, pfImpl, seq));
	}

	public double getKStarScoreLog10(int sequenceIndex, boolean useUB) {
		PFAbstract l = getPartitionFunction(1, getSequences(1).get(sequenceIndex));
		PFAbstract p = getPartitionFunction(0, getSequences(0).get(sequenceIndex));
		PFAbstract pl = getPartitionFunction(2, getSequences(2).get(sequenceIndex));
		return KSCalc.getKStarScoreLog10(l, p, pl, useUB);
	}
}
