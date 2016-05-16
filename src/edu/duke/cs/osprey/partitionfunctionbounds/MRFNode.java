/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.PositionConfSpaceSuper;
import edu.duke.cs.osprey.confspace.SuperRCTuple;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class MRFNode {

    //identifier of node, also accesses ematrix as posNum
    int posNum;
    //index into nodeList
    int index;
    //list of labels for the node; labels are RCs
    ArrayList<MRFLabel> labelList;
    int[] labels;
    
    int numLabels;
    //list of nodes that are neighbors of this node
    ArrayList<MRFNode> neighborList;
    public MRFNode(int posNum, ArrayList<Integer> unprunedRCs, int index) {
        this.posNum = posNum;
        this.index = index;
        //create label
        this.labelList = new ArrayList<>();
        this.labels = new int[unprunedRCs.size()];
        for (int i =0; i<unprunedRCs.size(); i++){
            int rc = unprunedRCs.get(i);
            MRFLabel label = new MRFLabel(rc);
            labelList.add(label);
            labels[i] = rc;
        }
        numLabels = labelList.size();
    }

    @Override
    public int hashCode() {
        int hash = 7;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final MRFNode other = (MRFNode) obj;
        if (this.posNum != other.posNum) {
            return false;
        }
        return true;
    }
    
    public int getNumLabels(){
        return this.numLabels;
    }

    
    
}
