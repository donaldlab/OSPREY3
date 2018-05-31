package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

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
	private MSSearchProblem search;//search problem
	private PruningMatrix pmat;//pruning matrix
	private HashMap<Integer, Integer> index2AbsolutePos;//maps from index to absolute position space
	//map from index to absolutePos. a* only sees index space. map back 
	//to absolute pos only when accessing the energy matrix or pruning matrix
	private HashMap<Integer, Integer> absolutePos2Index;//maps from absolute pos to index
	private RCTuple absoluteTuple;//tuple with positions converted to abs pos
	private Integer[] notAllowedPos;//unallowed positions
	private int[] totUndefinedPos;//undefined+notallowed
	private int totNumUndefined;
	private final int defined = -2;

	public MultiSequenceConfTree(MSSearchProblem search, EnergyMatrix emat, PruningMatrix pmat) {
		super(new FullAStarNode.Factory(search.getNumAssignedPos()), search, pmat);
		this.energyLBs = search.settings.energyLBs;
		this.search = search;
		this.emat = emat;
		this.pmat = pmat;
		this.index2AbsolutePos = new HashMap<Integer, Integer>();
		this.absolutePos2Index = new HashMap<Integer, Integer>();
		this.absoluteTuple = new RCTuple();
		init();
		setVerbose(false);
	}

	protected Integer[] getPosNums(boolean defined) {
		ArrayList<Integer> ans = search.getPosNums(defined);
		return ans.toArray(new Integer[ans.size()]);
	}

	protected void init() {
		totUndefinedPos = new int[search.confSpace.numPos];

		Integer[] allowedPos = getPosNums(true);
		numPos = allowedPos.length;

		//map from index to absolutePos. a* only sees index space. map back 
		//to absolute pos only when accessing the energy matrix or pruning matrix
		for(int i=0;i<numPos;++i) index2AbsolutePos.put(i, allowedPos[i]);
		notAllowedPos = getPosNums(false);
		for(int i=0;i<notAllowedPos.length;++i) index2AbsolutePos.put(i+numPos, notAllowedPos[i]);
		for(int key : index2AbsolutePos.keySet()) absolutePos2Index.put(index2AbsolutePos.get(key), key);

		assert(index2AbsolutePos.size()==search.confSpace.numPos);

		definedPos = new int[numPos];
		definedRCs = new int[numPos];
		undefinedPos = new int[numPos];
		childConf = new int[numPos];

		// see which RCs are unpruned and thus available for consideration
		// pack them into an efficient int matrix
		unprunedRCsAtPos = new int[search.confSpace.numPos][];
		for (int pos=0;pos<unprunedRCsAtPos.length;++pos) {//store all assigned and unassigned
			ArrayList<Integer> srcRCs = pmat.unprunedRCsAtPos(index2AbsolutePos.get(pos));
			int[] destRCs = new int[srcRCs.size()];
			for (int i=0; i<srcRCs.size(); i++) destRCs[i] = srcRCs.get(i);
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

	//convert to absolute positions
	private RCTuple setAbsolutePos(RCTuple other) {
		absoluteTuple.set(other);
		for(int i=0;i<absoluteTuple.size();++i) {
			absoluteTuple.pos.set(i, index2AbsolutePos.get(absoluteTuple.pos.get(i)));
		}
		return absoluteTuple;
	}
	
	protected double scoreNode(int childPos, int childRc, int[] conf) {		
		if(traditionalScore) {
			rcTuple.set(conf);
			//score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
			//plus contributions associated with each of the undefined res ("h-score")

			//g-score
			absoluteTuple = setAbsolutePos(rcTuple);
			double gScore = emat.getConstTerm() + emat.getInternalEnergy(absoluteTuple);//intra+pairwise
			
			boolean fullyAssigned = true;
			int len = numPos;
			for(int i=0;i<len;++i) {
				if(conf[i]==-1) fullyAssigned = false;
			}
			
			//h-score
			//first fill in totUndefined
			int k, l; totNumUndefined = 0;
			for (k=0; k<numUndefined; k++) totUndefinedPos[k] = undefinedPos[k];
			for (l=0; l<notAllowedPos.length; ++l) totUndefinedPos[k+l] = absolutePos2Index.get(notAllowedPos[l]);
			Arrays.fill(totUndefinedPos, k+l, totUndefinedPos.length, defined);
			totNumUndefined = k+l;
			
			double hScore = 0;
			for (int pos1 : totUndefinedPos) {
				//skip if defined in parent
				if(pos1 == defined) continue;

				//skip if defined in child
				if (pos1 == childPos) continue;

				//bound on contribution of this residue
				double posBestE = fullyAssigned && !energyLBs ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
				int[] rc1s = unprunedRCsAtPos[pos1];
				int n1 = rc1s.length;
				for (int j=0; j<n1; j++) {
					int rc1 = rc1s[j];
					// OPTIMIZATION: manually inlining this is noticeably slower
					// maybe it causes too much register pressure
					double undefE = getUndefinedRCEnergy(conf, fullyAssigned, pos1, rc1, j, childPos, childRc);
					posBestE = fullyAssigned && !energyLBs ? Math.max(posBestE, undefE) : Math.min(posBestE, undefE);
				}

				hScore += posBestE;
			}
			
			if(!energyLBs) {
				double minValue = -1.0*Double.MAX_VALUE;
				if(gScore > 0) gScore = Double.isInfinite(gScore) ? Double.MAX_VALUE : gScore;
				else gScore = Double.isInfinite(gScore) ? minValue : gScore;
				
				if(hScore > 0) hScore = Double.isInfinite(hScore) ? Double.MAX_VALUE : hScore;
				else hScore = Double.isInfinite(hScore) ? minValue : hScore;
				
				if(gScore==Double.MAX_VALUE && hScore==gScore) return gScore;
				else if(gScore==minValue && hScore==gScore) return gScore;
				else if(gScore==Double.MAX_VALUE || hScore==Double.MAX_VALUE) return Double.MAX_VALUE;
			}
			
			return gScore+hScore;

		} else {
			//other possibilities include MPLP, etc.
			//But I think these are better used as refinements
			//we may even want multiple-level refinement
			throw new UnsupportedOperationException("Advanced A* scoring methods not implemented yet!");
		}
	}

	private double getUndefinedRCEnergy(int[] conf, boolean fullyAssigned, int pos1, int rc1, int rc1i, int childPos, int childRc) {
		assertSplitPositions();
		//Provide a lower bound on what the given rc at the given level can contribute to the energy
		//assume partialConf and definedTuple

		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations (like preferring stack over heap) can make an impact

		// that said, let's copy some references to the stack =)
		EnergyMatrix emat = this.emat;
		int numDefined = this.numDefined;
		int numUndefined = this.totNumUndefined;
		int[] definedPos = this.definedPos;
		int[] definedRCs = this.definedRCs;
		int[] undefinedPos = this.totUndefinedPos;

		double rcContrib = emat.getOneBody(index2AbsolutePos.get(pos1), rc1);

		//for this kind of lower bound, we need to split up the energy into the defined-tuple energy
		//plus "contributions" for each undefined residue
		//so we'll say the "contribution" consists of any interactions that include that residue
		//but do not include higher-numbered undefined residues

		// first pass, defined residues
		for (int i=0; i<numDefined; i++) {
			int pos2 = definedPos[i];
			int rc2 = definedRCs[i];

			assert (pos2 != childPos);

			rcContrib += emat.getPairwise(index2AbsolutePos.get(pos1), rc1, index2AbsolutePos.get(pos2), rc2);
			//add higher-order terms that involve rc, rc2, and parts of partialConf
			//besides that only residues in definedTuple or levels below pos2
			//rcContrib += higherOrderContribLB(conf, pos1, rc1, pos2, rc2);
		}

		// if the child has a new definition, add that too
		if (childPos >= 0) {
			rcContrib += emat.getPairwise(index2AbsolutePos.get(pos1), rc1, index2AbsolutePos.get(childPos), childRc);
			//rcContrib += higherOrderContribLB(conf, pos1, rc1, childPos, childRc);
		}

		// second pass, undefined residues
		for (int i=0; i<numUndefined; i++) {
			int pos2 = undefinedPos[i];
			if (pos2 >= pos1) {
				break;
			}

			// skip if defined in child
			if (pos2 == childPos) {
				continue;
			}

			// min/max over all possible conformations
			double bestE = fullyAssigned && !energyLBs ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
			for (int rc2 : this.unprunedRCsAtPos[pos2]) {
				double pairwiseEnergy = emat.getPairwise(index2AbsolutePos.get(pos1), rc1, index2AbsolutePos.get(pos2), rc2);
				//pairwiseEnergy += higherOrderContribLB(conf, pos1, rc1, pos2, rc2);
				bestE = fullyAssigned && !energyLBs ? Math.max(bestE, pairwiseEnergy) : Math.min(bestE, pairwiseEnergy);
			}

			rcContrib += bestE;
		}

		return rcContrib;
	}

	protected double scoreNode(int[] conf) {
		if(traditionalScore) {
			rcTuple.set(conf);
			absoluteTuple = setAbsolutePos(rcTuple);

			//score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
			//plus contributions associated with each of the undefined res ("h-score")

			//"g-score"
			double gScore = emat.getConstTerm() + emat.getInternalEnergy(absoluteTuple);//intra+pairwise
			//defined to definåed energies

			//"h-score"
			boolean fullyAssigned = true;
			int len = numPos;
			for(int i=0;i<len;++i) {
				if(conf[i]==-1) fullyAssigned = false;
			}
			
			double hScore = 0;
			for(int pos=0; pos<search.confSpace.numPos;++pos) {
				if(rcTuple.pos.contains(pos)) continue;//already computed internal energy for defined pos
				double posBestE = fullyAssigned && !energyLBs ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
				for(int rc : unprunedRCsAtPos[pos]) {//pos here is undefined: either yet to be assigned or multi-sequence
					double undefE = getUndefinedRCEnergy(pos, rc, rcTuple, fullyAssigned);
					posBestE = fullyAssigned && !energyLBs ? Math.max(posBestE, undefE) : Math.min(posBestE, undefE);
				}

				hScore += posBestE;
			}
			return gScore + hScore;
		}

		else {
			//other possibilities include MPLP, etc.
			//But I think these are better used as refinements
			//we may even want multiple-level refinement
			throw new UnsupportedOperationException("Advanced A* scoring methods not implemented yet!");
		}
	}

	protected double getUndefinedRCEnergy(int pos1, int rc1, RCTuple definedTuple, boolean fullyAssigned) {
		//Provide a lower bound on what the given rc at the given level can contribute to the energy
		//assume partialConf and definedTuple

		//pos 1 is undefined: either yet to be assigned or multi-sequence;
		double rcContrib = emat.getOneBody(index2AbsolutePos.get(pos1), rc1);

		//for this kind of lower bound, we need to split up the energy into the defined-tuple energy
		//plus "contributions" for each undefined residue
		//so we'll say the "contribution" consists of any interactions that include that residue
		//but do not include higher-numbered undefined residues
		for(int pos2 = 0; pos2 < search.confSpace.numPos; pos2++){

			if(definedTuple.pos.contains(pos2) || pos2 < pos1) {//defined or lower numbered residues

				double posBestE = fullyAssigned && !energyLBs ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;//best pairwise energy

				for(int rc2 : unprunedRCsAtPos[pos2]) {

					double pairWiseE = emat.getPairwise(index2AbsolutePos.get(pos1), rc1, index2AbsolutePos.get(pos2), rc2);
					double higherOrderE = higherOrderContrib(pos1, rc1, pos2, rc2, definedTuple);
					pairWiseE += higherOrderE;

					posBestE = fullyAssigned && !energyLBs ? Math.max(posBestE, pairWiseE) : Math.min(posBestE, pairWiseE);
				}

				rcContrib += posBestE;
			}
		}

		return rcContrib;
	}

	protected double higherOrderContrib(int pos1, int rc1, int pos2, int rc2, RCTuple definedTuple) {
		//higher-order contribution for a given RC pair, when scoring a partial conf

		HigherTupleFinder<Double> htf = emat.getHigherOrderTerms(index2AbsolutePos.get(pos1), rc1, index2AbsolutePos.get(pos2), rc2);

		if(htf==null)
			return 0;//no higher-order interactions
		else {
			RCTuple curPair = new RCTuple(index2AbsolutePos.get(pos1), rc1, index2AbsolutePos.get(pos2), rc2);
			absoluteTuple = setAbsolutePos(definedTuple);
			return higherOrderContrib(htf, curPair, absoluteTuple);
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

		//double ans = scoreNode(conf);
		double ans = scoreNode(childPos, childRc, conf);
		return ans;
	}

}
