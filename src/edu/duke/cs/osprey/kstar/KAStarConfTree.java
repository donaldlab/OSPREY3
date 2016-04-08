package edu.duke.cs.osprey.kstar;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */

/*
 * todo
 * 1) get g-score wrt pan seq?
 */

@SuppressWarnings("serial")
public class KAStarConfTree extends ConfTree implements Serializable, ConfSearch {

	protected SearchProblem panSeqSP = null;
	protected HashMap<Integer, Integer> seq2PanSeq = null;
	protected HashMap<Integer, Integer> panSeq2Seq = null;
	protected boolean energyLB = true; // compute lb xor ub
	protected boolean usePrunedConfs = false;


	public KAStarConfTree(SearchProblem sp, SearchProblem panSeqSP, 
			ArrayList<Integer> panSeqPos, boolean usePrunedConfs) {
		super(sp);
		
		this.panSeqSP = panSeqSP;
		mapSeq2PanSeq(panSeqPos);
		energyLB = panSeqSP.contSCFlex ? true : false;
		this.usePrunedConfs = usePrunedConfs;
		
		init(sp, sp.pruneMat);
	}

    protected void init(SearchProblem sp, PruningMatrix pruneMat) {
    	unprunedRCsAtPos.clear();
        
        //see which RCs are unpruned and thus available for consideration
        for(int pos=0; pos<numPos; pos++){
            unprunedRCsAtPos.add(getRCsAtPos(pruneMat, pos));
        }
    }
	
    
	protected void mapSeq2PanSeq(ArrayList<Integer> panSeqPos) {
		
		seq2PanSeq = new HashMap<>(numPos);
		panSeq2Seq = new HashMap<>(numPos);
		
		for(int i = 0; i < numPos; ++i) {
			seq2PanSeq.put(i, panSeqPos.get(i));
			panSeq2Seq.put(panSeqPos.get(i), i);
		}
	}
	
	
	protected ArrayList<Integer> getPanSeqPos( RCTuple definedTuple ) {
		ArrayList<Integer> ans = new ArrayList<>(definedTuple.pos.size());
		for( int i : definedTuple.pos ) ans.add( seq2PanSeq.get(i) );
		return ans;
	}
	

	protected ArrayList<Integer> getUndefinedPos( RCTuple definedTuple ) {
		ArrayList<Integer> ans = new ArrayList<>(panSeqSP.confSpace.numPos);
		
		for(int level = 0; level < panSeqSP.confSpace.numPos; ++level) ans.add(level);

		ArrayList<Integer> definedPos = getPanSeqPos(definedTuple);
		
		// remove defined positions
		ans.removeAll(definedPos);

		return ans;
	}


	protected ArrayList<Integer> allowedRCsAtLevel(int level, int[] partialConf, ArrayList<Integer> undefinedPos) {
		// What RCs are allowed at the specified level (i.e., position num) in the given partial conf?
		ArrayList<Integer> allowedRCs;

		if(undefinedPos.contains(level)) { // level is undefined
			allowedRCs = getRCsAtPos(panSeqSP.pruneMat, level);
		}

		else {
			allowedRCs = new ArrayList<>();
			allowedRCs.add(partialConf[panSeq2Seq.get(level)]);
		}

		return allowedRCs;
	}


	protected double RCContribution(int level, int rc, RCTuple definedTuple, int[] partialConf, 
			ArrayList<Integer> undefinedPos) {
		//Provide a lower bound on what the given rc at the given level can contribute to the energy
		//assume partialConf and definedTuple

		double rcContrib = panSeqSP.getEnergyMatrix().getOneBody(level, rc);

		//for this kind of lower bound, we need to split up the energy into the defined-tuple energy
		//plus "contributions" for each undefined residue
		//so we'll say the "contribution" consists of any interactions that include that residue
		//but do not include higher-numbered undefined residues
		for(int level2 = 0; level2 < panSeqSP.confSpace.numPos; level2++){

			if(!undefinedPos.contains(level2) || level2 < level){//defined or lower numbered residues

				double levelBestE = energyLB ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;//best pairwise energy

				ArrayList<Integer> allowedRCs = allowedRCsAtLevel(level2, partialConf, undefinedPos);

				for( int rc2 : allowedRCs ) {
					
					if(!panSeqSP.pruneMat.getPairwise(level, rc, level2, rc2)) {

						double interactionE = panSeqSP.getEnergyMatrix().getPairwise(level, rc, level2, rc2);

						double higherOrderE = higherOrderContrib(level, rc, level2, rc2, partialConf, undefinedPos);
						//add higher-order terms that involve rc, rc2, and parts of partialConf

						interactionE += higherOrderE;

						//besides that only residues in definedTuple or levels below level2
						levelBestE = energyLB ? Math.min(levelBestE, interactionE) : Math.max(levelBestE, interactionE);
					}
				}

				rcContrib += levelBestE;
			}
		}

		return rcContrib;
	}


	protected double higherOrderContrib(int pos1, int rc1, int pos2, int rc2, 
			int[] partialConf, ArrayList<Integer> undefinedPos) {
		//higher-order contribution for a given RC pair, when scoring a partial conf

		HigherTupleFinder<Double> htf = panSeqSP.emat.getHigherOrderTerms(pos1, rc1, pos2, rc2);

		if(htf==null)
			return 0;//no higher-order interactions
		else {
			RCTuple curPair = new RCTuple(pos1, rc1, pos2, rc2);
			return higherOrderContrib(htf, curPair, partialConf, undefinedPos);
		}
	}


	@SuppressWarnings("unchecked")
	double higherOrderContrib(HigherTupleFinder<Double> htf, RCTuple startingTuple, int[] partialConf, ArrayList<Integer> undefinedPos) {
		//recursive function to get bound on higher-than-pairwise terms
		//this is the contribution to the bound due to higher-order interactions
		//of the RC tuple startingTuple (corresponding to htf)

		double contrib = 0;

		//to avoid double-counting, we are just counting interactions of starting tuple
		//with residues before the "earliest" one (startingLevel) in startingTuple
		//"earliest" means lowest-numbered, except non-mutating res come before mutating
		int startingLevel = startingTuple.pos.get( startingTuple.pos.size()-1 );

		for(int iPos : htf.getInteractingPos()){//position has higher-order interaction with tup
			if(posComesBefore(iPos,startingLevel,undefinedPos)) {//interaction in right order
				//(want to avoid double-counting)

				double levelBestE = energyLB ? Double.POSITIVE_INFINITY : 
					Double.NEGATIVE_INFINITY;//best value of contribution from tup-iPos interaction

				ArrayList<Integer> allowedRCs = allowedRCsAtLevel(iPos, partialConf, undefinedPos);

				for( int rc : allowedRCs ){

					RCTuple augTuple = startingTuple.addRC(iPos, rc);

					if( !panSeqSP.pruneMat.isPruned(augTuple) ){

						double interactionE = htf.getInteraction(iPos, rc);

						//see if need to go up to highers order again...
						@SuppressWarnings("rawtypes")
						HigherTupleFinder htf2 = htf.getHigherInteractions(iPos, rc);
						if(htf2!=null){
							interactionE += higherOrderContrib(htf2, augTuple, partialConf, undefinedPos);
						}

						levelBestE = energyLB ? Math.min(levelBestE, interactionE) : Math.max(levelBestE, interactionE);
					}
				}

				contrib += levelBestE;//add up contributions from different interacting positions iPos
			}
		}

		return contrib;
	}


	protected boolean posComesBefore(int pos1, int pos2, ArrayList<Integer> undefinedPos){
		//for purposes of contributions to traditional conf score, 
		//we go through defined and then through undefined positions (in partialConf);
		//within each of these groups we go in order of position number
		if(!undefinedPos.contains(pos2)){//pos2 defined
			return (pos1<pos2 && !undefinedPos.contains(pos1));//pos1 must be defined to come before pos2
		}
		else//pos1 comes before pos2 if it's defined, or if pos1<pos2
			return (pos1<pos2 || !undefinedPos.contains(pos1));
	}


	protected double scoreConf(int[] partialConf){
		
		if(traditionalScore) {
			RCTuple definedTuple = new RCTuple(partialConf);

			//"g-score"
			double score = emat.getConstTerm() + emat.getInternalEnergy(definedTuple);

			//"h-score"
			//score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
			//plus contributions associated with each of the undefined res ("h-score")
			ArrayList<Integer> undefinedPos = getUndefinedPos(definedTuple);

			for(int level : undefinedPos) {

				double bestInteractionE = energyLB ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;

				ArrayList<Integer> rotList = getRCsAtPos(panSeqSP.pruneMat, level);
				for(int rc : rotList) {
					double rcContribution = RCContribution(level, rc, definedTuple, partialConf, undefinedPos);

					bestInteractionE = energyLB ? Math.min(bestInteractionE, rcContribution) : 
						Math.max(bestInteractionE, rcContribution);
				}

				score += bestInteractionE;
			}

			return score;
		}

		else {
			//other possibilities include MPLP, etc.
			//But I think these are better used as refinements
			//we may even want multiple-level refinement
			throw new RuntimeException("Advanced A* scoring methods not implemented yet!");
		}
	}
	
	
	public double confBound(int[] partialConf) {
		return scoreConf(partialConf);
	}
	
	
	ArrayList<Integer> getRCsAtPos(PruningMatrix pruneMat, int pos) {
		return usePrunedConfs ? pruneMat.prunedRCsAtPos(pos) : pruneMat.unprunedRCsAtPos(pos);
	}
}
