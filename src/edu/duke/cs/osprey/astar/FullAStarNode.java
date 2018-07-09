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

package edu.duke.cs.osprey.astar;

import java.io.Serializable;
import java.util.Arrays;

public class FullAStarNode implements AStarNode, Serializable {
	
	private static final long serialVersionUID = -537132381411057989L;

	public static class Factory implements AStarNode.Factory<FullAStarNode>, Serializable {

		private static final long serialVersionUID = -6909420740336320965L;
		private int numPos;
		
		public Factory(int numPos) {
			this.numPos = numPos;
		}
		
		@Override
		public FullAStarNode makeRoot() {
			int conf[] = new int[numPos];
			Arrays.fill(conf, -1);
			return new FullAStarNode(conf);
		}
		
		@Override
		public FullAStarNode make(FullAStarNode parent, int pos, int rc) {
			
			// explicitly instantiate the conformation
            int[] conf = parent.getNodeAssignments().clone();
            conf[pos] = rc;
            
			return new FullAStarNode(conf);
		}
	}
	
    private int nodeAssignments[];//assignments (e.g. partial conformation) for node
    
    private double score;//score (probably a lower bound on the energy)
    private double gscore;
    private double hscore;
    
    //indicates the score needs to be refined (e.g. with EPIC continuous terms)
    //always false in simpler versions of A*
    boolean scoreNeedsRefinement;

    
    //These are used in COMETS
    public double UB = Double.POSITIVE_INFINITY;//upper bound
    public int UBConf[] = null;//can have an upper bound on GMEC energy for this node's conf space
    //(and thus on the overall GMEC energy)
    //that is the energy of the conf denoted by UBConf
    
    
    

    public FullAStarNode(int[] nodeAssignments) {
        this.nodeAssignments = nodeAssignments;
        this.score = Double.NaN;
        this.gscore = Double.NaN;
        this.hscore = Double.NaN;
        this.scoreNeedsRefinement = false;
    }
    
    public FullAStarNode(FullAStarNode fan){
        //copy constructor.  Shallow copy
        UB = fan.UB;
        UBConf = fan.UBConf;
        gscore = fan.gscore;
        hscore = fan.hscore;
        nodeAssignments = fan.nodeAssignments;
        score = fan.score;
        scoreNeedsRefinement = fan.scoreNeedsRefinement;
    }

    @Override
    public int compareTo(AStarNode other) {
        return Double.valueOf(score).compareTo(other.getScore());
    }

    @Override
    public int[] getNodeAssignments() {
        return nodeAssignments;
    }
    
    @Override
    public void setScore(double score) {
        this.score = score;
    }
    
    @Override
    public double getScore() {
        return score;
    }
    
    public double getGScore() {
    	return gscore;
    }
    public void setGScore(double val) {
    	gscore = val;
    }
    
    public double getHScore() {
    	return hscore;
    }
    public void setHScore(double val) {
    	hscore = val;
    }
    
    @Override
    public int getLevel() {
        int level = 0;
        for (int a : nodeAssignments) {
            if (a >= 0) {
                level++;
            }
        }
        return level;
    }
    
    @Override
    public boolean isFullyDefined() {
        //Assuming assignments greater than 0 denote fully defined positions,
        //determine if this node is fully defined or not
        for(int a : nodeAssignments){
            if(a<0)
                return false;
        }
        return true;
    }

	@Override
	public boolean scoreNeedsRefinement() {
		return scoreNeedsRefinement;
	}
	
	@Override
	public void setScoreNeedsRefinement(boolean val) {
		scoreNeedsRefinement = val;
	}
}
