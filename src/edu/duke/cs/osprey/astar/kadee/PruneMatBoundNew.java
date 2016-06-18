/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.kadee;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author hmn5
 */
public class PruneMatBoundNew extends PruningMatrix {

    List<Boolean> isUnassigned;
    PruningMatrix parent;

    public PruneMatBoundNew(List<Boolean> isUnassigned, PruningMatrix parent) {
        this.isUnassigned = isUnassigned;
        this.parent = parent;
    }

    @Override
    public Boolean getOneBody(int res, int index) {
        if (isUnassigned.get(res)) {
            return false;
        } else {
            return parent.getOneBody(res, index);
        }
    }

    @Override
    public Boolean getPairwise(int res1, int index1, int res2, int index2) {
        if (isUnassigned.get(res1) && isUnassigned.get(res2)) {
            return false;
        }
        return parent.getPairwise(res1, index1, res2, index2);
    }

    @Override
    public int numRCsAtPos(int pos) {
        return parent.numRCsAtPos(pos);
    }

    @Override
    public ArrayList<Integer> unprunedRCsAtPos(int pos) {
        //which RCs at the given position are unpruned?
        //Return index of the RCs within the position
        ArrayList<Integer> ans = new ArrayList<>();
        int numRCs = numRCsAtPos(pos);

        for (int index = 0; index < numRCs; index++) {
            if (!getOneBody(pos, index)) {
                ans.add(index);
            }
        }

        return ans;
    }
    
    @Override
    public boolean isPruned(RCTuple tup) {
        //can be prune per se, or check if some singles in it are pruned, or pairs, etc.
        for (int indexInTup = 0; indexInTup < tup.pos.size(); indexInTup++) {
            int pos1 = tup.pos.get(indexInTup);
            int rc1 = tup.RCs.get(indexInTup);

            if (getOneBody(pos1, rc1))//rc1 pruned by itself
            {
                return true;
            }

            //check if any pairs pruned
            for (int index2 = 0; index2 < indexInTup; index2++) {
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);

                if (getPairwise(pos1, rc1, pos2, rc2)) {
                    return true;
                }
            }
        }
        return false;
    }

}
