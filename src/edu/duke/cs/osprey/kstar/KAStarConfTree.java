package edu.duke.cs.osprey.kstar;

import java.io.Serializable;
import java.util.ArrayList;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

@SuppressWarnings("serial")
public class KAStarConfTree extends ConfTree implements Serializable {

	protected EnergyMatrix allSeqEmat = null;
	protected ArrayList<ArrayList<Integer>> allSeqUnprunedRCsAtPos = new ArrayList<>();
	// compute lb xor ub
	protected boolean computeLB = true;

	public KAStarConfTree(SearchProblem sp, SearchProblem allSeqSP) {
		super(sp);
		
		allSeqInit(allSeqSP, allSeqSP.pruneMat);
	}


	public KAStarConfTree(SearchProblem sp, PruningMatrix pruneMat, boolean useEPIC, SearchProblem allSeqSP){
		super(sp,pruneMat,useEPIC);
		
		allSeqInit(allSeqSP, allSeqSP.pruneMat);
	}


	private void allSeqInit(SearchProblem allSeqSP, PruningMatrix allSeqpruneMat) {

		computeLB = allSeqSP.contSCFlex ? true : false;
		int allSeqNumPos = allSeqSP.confSpace.numPos;
		
		//see which RCs are unpruned and thus available for consideration
		for(int pos=0; pos<allSeqNumPos; pos++){
			allSeqUnprunedRCsAtPos.add( allSeqpruneMat.unprunedRCsAtPos(pos) );
		}

		//get the appropriate energy matrix to use in this A* search
		if(allSeqSP.useTupExpForSearch)
			allSeqEmat = allSeqSP.tupExpEMat;

		else
			allSeqEmat = allSeqSP.emat;
	}


	protected double scoreConf(int[] partialConf){
		if(traditionalScore){
			RCTuple definedTuple = new RCTuple(partialConf);

			//"g-score"
			double score = emat.getConstTerm() + emat.getInternalEnergy( definedTuple );

			//score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
			//plus contributions associated with each of the undefined res ("h-score")
			for(int level=0; level<numPos; level++){
				if(partialConf[level]<0){//level not fully defined

					double resContribLB = Double.POSITIVE_INFINITY;//lower bound on contribution of this residue
					//resContribLB will be the minimum_{rc} of the lower bound assuming rc assigned to this level

					for ( int rc : unprunedRCsAtPos.get(level) ) {
						resContribLB = Math.min(resContribLB, RCContributionLB(level,rc,definedTuple,partialConf));
					}

					score += resContribLB;
				}
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
}
