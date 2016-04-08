package edu.duke.cs.osprey.kstar;

import java.io.Serializable;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.PruningMatrix;

@SuppressWarnings("serial")
public class PStarConfTree extends ConfTree implements Serializable, ConfSearch{

	public PStarConfTree(SearchProblem sp) {
		super(sp);
		
		init(sp, sp.pruneMat);
	}
	
	protected void init(SearchProblem sp, PruningMatrix pruneMat) {
    	unprunedRCsAtPos.clear();
        
        //see which RCs are unpruned and thus available for consideration
        for(int pos=0; pos<numPos; pos++){
            unprunedRCsAtPos.add(pruneMat.prunedRCsAtPos(pos));
        }
    }
	
	public double confBound(int[] conf) {
		return scoreConf(conf);
	}
	
}
