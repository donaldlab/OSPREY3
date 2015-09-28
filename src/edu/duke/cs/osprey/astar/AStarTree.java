/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.confspace.ConfSearch;
import java.util.ArrayList;
import java.util.PriorityQueue;

/**
 *
 * @author mhall44
 */

public abstract class AStarTree implements ConfSearch {
    //This replaces MSAStar with something more generic
    //The goal is that subclasses of this, differing only in the abstract methods,
    //can cover the A* variations we're considering:
    //dynamic ordering, different node scoring heuristics, COMETS, super-rotamers, etc.
        
    private PriorityQueue<AStarNode> pq = null;
        
    //AStarNode can be lightweight: just int[], score, and flag for if score needs refinement
    //the meanings are assigned by subclasses of this class, which define things like scoring
    //and thus what the int[] means
    //Methods like COMETS can of course subclass AStarNode to include more information in the node
            
    public int numExpanded = 0;//counting number of nodes expanded
    public int numPruned = 0;//counting number of nodes pruned
    
    @Override
    public int[] nextConf(){
        //return best conformation remaining in tree
        
        if(pq==null){//need to initialize tree (indicates haven't enumerated anything from this tree yet)
            initQueue(rootNode());
        }
        
        AStarNode curNode;
        
        while(true) {//keep going until we either find the optimal solution, or find the tree is empty
            curNode = pq.poll();
            
            if(curNode==null){
                System.out.println("A* tree empty...returning empty signal");
                return null;//signal for empty tree
            }
            
            if(canPruneNode(curNode))//like, too many AA changes
                numPruned++;
            else {
                
                while(curNode.scoreNeedsRefinement){
                    refineScore(curNode);
                    if(curNode.score!=Double.POSITIVE_INFINITY)//remove node if refinement showed it's impossible
                        pq.add(curNode);
                    
                    curNode = pq.poll();
                    if(curNode==null){
                        System.out.println("A* tree empty...returning empty signal");
                        return null;//signal for empty tree
                    }
                }
                
                if(isFullyAssigned(curNode)){
                    return outputNode(curNode);
                }

                //expand
                ArrayList<AStarNode> children = getChildren(curNode);
                //note: in a method like COMETS that refines nodes,
                //expandNode may return a singleton list consisting of curNode with improved bound

                numExpanded++;
                
                for(AStarNode child : children)
                    pq.add(child);
            }
        }
        
    }
    
    
    public void initQueue(AStarNode node){
        pq = new PriorityQueue<>();
        pq.add(node);
    }
    
    
    //methods with default implementations that may need to be overridden:
    
    public boolean canPruneNode(AStarNode node){
        //By default we don't have node pruning
        //but subclasses may allow this
        return false;
    }
    
    
    public int[] outputNode(AStarNode node){
        //by default, the output of the A* tree will be simply the node assignments for the optimal node
        //but we may sometimes want to process it in some way
        System.out.println("A* returning conf.  "+pq.size()+" nodes in A* tree.  Score: "+node.score);
        return node.nodeAssignments;
    }
    
    
    void refineScore(AStarNode node){//e.g. add the EPIC contribution
        //In trees without score refinement, we should 
        //not be calling this method
        throw new UnsupportedOperationException("ERROR: Score refinement not supported"
                + " in this type of A* tree");
    }
  
    //abstract methods:
    
    
    //getChildren and rootNode create nodes with a score
    //the score can be a quick, possibly loose lower bound for now; 
    //mark scoreNeedRefinement if we'll want to refine any nodes that get to be head node
    
    public abstract ArrayList<AStarNode> getChildren(AStarNode curNode);
        //Get children for a node
        //this can be either static or dynamic ordering, depending on implementation
    
    public abstract AStarNode rootNode();
    
    
    public abstract boolean isFullyAssigned(AStarNode node);//is the node fully assigned (i.e., returnable?)

    
    
    public PriorityQueue<AStarNode> getQueue() {
        //direct access to queue.  Use with caution.  
        return pq;
    }
    
    
    
        
}
