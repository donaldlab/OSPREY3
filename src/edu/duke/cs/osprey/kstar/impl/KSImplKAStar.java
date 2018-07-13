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

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import edu.duke.cs.osprey.kstar.KSAllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KAStarNode;
import edu.duke.cs.osprey.kstar.KAStarTree;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

public class KSImplKAStar extends KSAbstract {

	public static boolean useTightBounds = true;
	public static String nodeExpansionMethod = "parallel1";

	public KSImplKAStar(KSConfigFileParser cfp) {
		super(cfp);

		if(useVoxelG){
			throw new RuntimeException("ERROR: Rigid-rotamer upper bound used by KSImplKAStar "
					+ "not compatible with continuous entropy (useVoxelG");
		}
	}

	@Override
	public void init( HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs ) {

		this.strand2AllowedSeqs = strand2AllowedSeqs;
                checkAPPP();

		printSequences();

		createOutputDir();

		createEmatDir();

		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, false));
		createEmats(contSCFlexVals);
	}


	@Override
	public String getKSMethod() {
		return "kastar";
	}

	protected BigInteger countProcessedConfs() {
		BigInteger ans = BigInteger.ZERO;

		for(PFAbstract pf : name2PF.values()) { 
			if(pf.isFullyDefined() && pf.getImpl().equalsIgnoreCase(PFAbstract.getCFGImpl()))
				ans = ans.add(pf.getNumProcessed());
		}

		return ans;
	}


	@Override
	public void run() {

		setStartTime(System.currentTimeMillis());

		int numOutput = 0;

		if(useTightBounds)
			numOutput = runTB();

		else
			numOutput = runLB();

		System.out.println("\nseqsCompleted: " + KAStarNode.getNumLeavesCompleted() 
		+ " seqsCreated: " + KAStarNode.getNumLeavesCreated()
		+ " seqsPossible: " + strand2AllowedSeqs.get(2).getNumSeqs()
		+ " seqsInOutput: " + numOutput);

		System.out.println("K* leaf conformations processed: " + countProcessedConfs());
		System.out.println("K* running time: " + (System.currentTimeMillis()-getStartTime())/1000 + " seconds\n");

		abortPFs();
	}


	private int runTB() {

		int completed = 0;

		// compute wt sequence for reference
		wtKSCalc = computeWTCalc();
		//setBestCalc(wtKSCalc);

		if( strand2AllowedSeqs.get(2).getNumSeqs() <= 1 )
			return completed; // wt is sequence[0]

		// initialize KUStar tree
		KAStarTree tree = new KAStarTree(this, strand2AllowedSeqs, wtKSCalc);

		// add root node
		tree.add( new KAStarNode(null, null, true) );

		int target = cfp.params.getInt("KStarNumSeqs", 5);

		for( KAStarNode best = tree.poll(); best != null && completed < target; 
				best = tree.poll() ) {

			if( best.isFullyProcessed() ) {

				best.checkConsistency(best);

				//if(passesInterMutationPruning(best.lb))
				best.lb.printSummary( getOputputFilePath(), getStartTime(), getNumSeqsCreated(0), getNumSeqsCompleted(0) );

				completed++;

				continue;
			}

			ArrayList<KAStarNode> children = best.expand();
			tree.add(children);
		}

		return completed + 1; // + 1 to count wild type
	}


	private int runLB() {
		// run until the lower bound of the next completed sequence is greater 
		// than the upper bound of any previously completed sequence

		int completed = 0;

		// compute wt sequence for reference
		wtKSCalc = computeWTCalc();

		if( strand2AllowedSeqs.get(2).getNumSeqs() <= 1 )
			return completed; // wt is sequence[0]

		// initialize KUStar tree
		KAStarTree tree = new KAStarTree(this, strand2AllowedSeqs, wtKSCalc);

		// add root node
		tree.add( new KAStarNode(null, null, true) );

		double gUB = Double.NEGATIVE_INFINITY;

		for( KAStarNode best = tree.poll(); best != null; best = tree.poll() ) {

			if( best.isFullyProcessed() ) {

				double uB = best.getUBScore();

				if( completed++ == 0 ) gUB = uB;

				if( best.getLBScore() > gUB && gUB > Double.NEGATIVE_INFINITY ) 
					break;
				else
					best.lb.printSummary( getOputputFilePath(), getStartTime(), getNumSeqsCreated(0), getNumSeqsCompleted(0) );

				if( uB > gUB ) gUB = uB;

				continue;
			}

			ArrayList<KAStarNode> children = best.expand();
			tree.add(children);
		}

		return completed + 1; // +1 to count wild type
	}
}
