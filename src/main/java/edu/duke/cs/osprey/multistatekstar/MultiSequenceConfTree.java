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

package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import java.util.ArrayList;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.FullAStarNode;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 * Multi-Sequence conformation tree used to enumerate energies for partial
 * conformations.
 * 
 */

@SuppressWarnings("serial")
public class MultiSequenceConfTree extends ConfTree<FullAStarNode> {

	public boolean energyLBs;//compute either energy lower bound or upper bound
	MSSearchProblem search;//search problem
	PruningMatrix pmat;//pruning matrix
	Integer[] allowedPos;//largest set of positions allowed by the (partial) sequence

	public MultiSequenceConfTree(MSSearchProblem search, EnergyMatrix emat, PruningMatrix pmat) {
		super(new FullAStarNode.Factory(search.getNumAssignedPos()), search, pmat);
		this.energyLBs = search.settings.energyLBs;
		this.search = search;
		this.emat = emat;
		this.pmat = pmat;
		allowedPos = getPosNums(true);
		init();
	}

	protected Integer[] getPosNums(boolean defined) {
		ArrayList<Integer> ans = search.getPosNums(defined);
		return ans.toArray(new Integer[ans.size()]);
	}

	protected void init() {
		numPos = allowedPos.length;

		definedPos = new int[numPos];
		definedRCs = new int[numPos];
		undefinedPos = new int[numPos];
		childConf = new int[numPos];

		// see which RCs are unpruned and thus available for consideration
		// pack them into an efficient int matrix
		unprunedRCsAtPos = new int[search.confSpace.numPos][];
		for (int pos=0;pos<unprunedRCsAtPos.length;++pos) {//store all assigned and unassigned
			ArrayList<Integer> srcRCs = pmat.unprunedRCsAtPos(pos);
			int[] destRCs = new int[srcRCs.size()];
			for (int i=0; i<srcRCs.size(); i++) {
				destRCs[i] = srcRCs.get(i);
			}
			unprunedRCsAtPos[pos] = destRCs;
		}

		//get the appropriate energy matrix to use in this A* search
		if(search.useTupExpForSearch)
			emat = search.tupExpEMat;
		else {
			emat = search.emat;

			if(search.useEPIC){//include EPIC in the search
				useRefinement = true;
				epicMat = search.epicMat;
				//confSpace = search.confSpace;
				minPartialConfs = search.epicSettings.minPartialConfs;
			}
		}
	}

	@Override
	public BigInteger getNumConformations() {
		BigInteger num = BigInteger.valueOf(1);
		for (int pos : allowedPos) {
			num = num.multiply(BigInteger.valueOf(unprunedRCsAtPos[pos].length));
		}
		return num;
	}

	protected void splitPositions(FullAStarNode node) {

		// make sure we're not split already
		assert (numDefined == 0 && numUndefined == 0);

		int[] conf = node.getNodeAssignments();

		// split conformation into defined and undefined residues
		numDefined = 0;
		numUndefined = 0;
		for (int pos : allowedPos) {
			int rc = conf[pos];
			if (rc >= 0) {
				definedPos[numDefined] = pos;
				definedRCs[numDefined] = rc;
				numDefined++;
			} else {
				undefinedPos[numUndefined] = pos;
				numUndefined++;
			}
		}

		assert (numDefined + numUndefined == numPos);
	}

	protected double scoreNode(int[] partialConf) {
		if(traditionalScore) {
			rcTuple.set(partialConf);
			//"g-score"
			double score = emat.getConstTerm() + emat.getInternalEnergy(rcTuple);//intra+pairwise

			//"h-score"
			//score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
			//plus contributions associated with each of the undefined res ("h-score")

			for(int pos=0; pos<search.confSpace.numPos;++pos) {
				if(rcTuple.pos.contains(pos)) continue;//skip positions assigned in rc tuple

				double bestE = energyLBs ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;

				for(int rc : unprunedRCsAtPos[pos]) {
					double rcContrib = RCContribution(pos, rc, rcTuple);
					bestE = energyLBs ? Math.min(bestE, rcContrib) : Math.max(bestE, rcContrib);
				}
				score += bestE;
			}
			return score;
		}

		else {
			//other possibilities include MPLP, etc.
			//But I think these are better used as refinements
			//we may even want multiple-level refinement
			throw new UnsupportedOperationException("Advanced A* scoring methods not implemented yet!");
		}
	}

	protected double RCContribution(int pos1, int rc1, RCTuple definedTuple) {
		//Provide a lower bound on what the given rc at the given level can contribute to the energy
		//assume partialConf and definedTuple

		double rcContrib = emat.getOneBody(pos1, rc1);

		//for this kind of lower bound, we need to split up the energy into the defined-tuple energy
		//plus "contributions" for each undefined residue
		//so we'll say the "contribution" consists of any interactions that include that residue
		//but do not include higher-numbered undefined residues
		for(int pos2 = 0; pos2 < search.confSpace.numPos; pos2++){

			if(definedTuple.pos.contains(pos2) || pos2 < pos1) {//defined or lower numbered residues

				double posBestE = energyLBs ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;//best pairwise energy

				for(int rc2 : unprunedRCsAtPos[pos2]) {

					double interactionE = emat.getPairwise(pos1, rc1, pos2, rc2);
					double higherOrderE = higherOrderContrib(pos1, rc1, pos2, rc2, definedTuple);
					interactionE += higherOrderE;

					posBestE = energyLBs ? Math.min(posBestE, interactionE) : Math.max(posBestE, interactionE);
				}

				rcContrib += posBestE;
			}
		}

		return rcContrib;
	}

	protected double higherOrderContrib(int pos1, int rc1, int pos2, int rc2, RCTuple definedTuple) {
		//higher-order contribution for a given RC pair, when scoring a partial conf

		HigherTupleFinder<Double> htf = emat.getHigherOrderTerms(pos1, rc1, pos2, rc2);

		if(htf==null)
			return 0;//no higher-order interactions
		else {
			RCTuple curPair = new RCTuple(pos1, rc1, pos2, rc2);
			return higherOrderContrib(htf, curPair, definedTuple);
		}
	}

	@SuppressWarnings("unchecked")
	double higherOrderContrib(HigherTupleFinder<Double> htf, RCTuple startingTuple,
			RCTuple definedTuple) {
		//recursive function to get bound on higher-than-pairwise terms
		//this is the contribution to the bound due to higher-order interactions
		//of the RC tuple startingTuple (corresponding to htf)

		double contrib = 0;

		//to avoid double-counting, we are just counting interactions of starting tuple
		//with residues before the "earliest" one (startingLevel) in startingTuple
		//"earliest" means lowest-numbered, except non-mutating res come before mutating
		int startingLevel = startingTuple.pos.get( startingTuple.pos.size()-1 );

		for(int iPos : htf.getInteractingPos()){//position has higher-order interaction with tup
			if(posComesBefore(iPos,startingLevel,definedTuple)) {//interaction in right order
				//(want to avoid double-counting)

				double posBestE = energyLBs ? Double.POSITIVE_INFINITY : 
					Double.NEGATIVE_INFINITY;//best value of contribution from tup-iPos interaction

				for(int rc : unprunedRCsAtPos[iPos]) {

					RCTuple augTuple = startingTuple.addRC(iPos, rc);

					double interactionE = htf.getInteraction(iPos, rc);

					//see if need to go up to highers order again...
					@SuppressWarnings("rawtypes")
					HigherTupleFinder htf2 = htf.getHigherInteractions(iPos, rc);
					if(htf2!=null){
						interactionE += higherOrderContrib(htf2, augTuple, definedTuple);
					}

					posBestE = energyLBs ? Math.min(posBestE, interactionE) : Math.max(posBestE, interactionE);
				}

				contrib += posBestE;//add up contributions from different interacting positions iPos
			}
		}

		return contrib;
	}

	protected boolean posComesBefore(int pos1, int pos2, RCTuple definedTuple){
		//for purposes of contributions to traditional conf score, 
		//we go through defined and then through undefined positions (in partialConf);
		//within each of these groups we go in order of position number
		if(definedTuple.pos.contains(pos2)){//pos2 defined
			return (pos1<pos2 && definedTuple.pos.contains(pos1));//pos1 must be defined to come before pos2
		}
		else//pos1 comes before pos2 if it's defined, or if pos1<pos2
			return (pos1<pos2 || definedTuple.pos.contains(pos1));
	}

	protected double scoreConfDifferential(FullAStarNode parentNode, int childPos, int childRc) {
		assertSplitPositions();

		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations can make an impact

		// get the full conf, start with the parent first
		int[] conf = parentNode.getNodeAssignments();

		// but if this is actually a child node, switch to the child conf
		if (childPos >= 0) {

			// parent shouldn't be assigned here
			assert (conf[childPos] < 0);

			// make the child conf
			System.arraycopy(conf, 0, childConf, 0, numPos);
			childConf[childPos] = childRc;
			conf = childConf;
		}

		double ans = scoreNode(conf);
		return ans;
	}

}
