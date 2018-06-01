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

package edu.duke.cs.osprey.astar.ewakstar;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.FullAStarNode;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.util.ArrayList;
import java.util.PriorityQueue;

public class EWAKStarNode extends FullAStarNode{

    PruningMatrix pruneMat;//pruning matrix for the PL state

    //These things are only needed (and defined) for fully defined sequences
    ConfTree<FullAStarNode> stateTree = null;//tree searching conf space for state
    double stateUB;//upper bounds on GMEC for PL state
    
    public EWAKStarNode(int[] nodeAssignments, PruningMatrix pruneMat) {
        super(nodeAssignments);//score not assigned yet, and doesn't need refinement
        this.pruneMat = pruneMat;
    }

    public EWAKStarNode(EWAKStarNode en){
        //copy constructor.  COMETS-specific things are big so let's not copy them unnecessarily
        super(en);
        pruneMat = en.pruneMat;
        stateTree = en.stateTree;
        stateUB = en.stateUB;
    }

    void expandConfTree(){

        //given a stateTree, look at the current queue and then get all of the children from the best node and
        //expand that node - add those children to the "expansion" - your tree is now expanded a level.
        PriorityQueue<FullAStarNode> expansion = stateTree.getQueue();

        FullAStarNode bestNode = expansion.poll();
        ArrayList<FullAStarNode> children = stateTree.getChildren(bestNode);

        for(FullAStarNode child : children){
            child.UB = Double.NaN;//update UB later if needed
            child.UBConf = bestNode.UBConf;//store this as an initial value
            expansion.add(child);
        }
    }

}

