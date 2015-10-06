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
    int nodeNum;
    //list of labels for the node; labels are super-RCs
    ArrayList<MRFLabel> labelList;
    //list of nodes that are neighbors of this node
    ArrayList<MRFNode> neighborList;
    public MRFNode(int posNum, ArrayList<Integer> unprunedSuperRCs) {
        this.nodeNum = posNum;

        //create label
        this.labelList = new ArrayList<>();
        for (int superRC : unprunedSuperRCs){
            MRFLabel label = new MRFLabel(superRC);
            labelList.add(label);
        }
    }
    
}
