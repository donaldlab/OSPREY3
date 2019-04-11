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

package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.tools.HashCalculator;

import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class RCTuple implements Serializable {
    
	private static final long serialVersionUID = 3773470316855163174L;
	
	//a tuple of RCs
    public ArrayList<Integer> pos;//which flexible positions
    public ArrayList<Integer> RCs;//the RCs themselves (residue-specific numbering, as in the TupleMatrices)

    
    public RCTuple(){
        //empty pos, RCs (basically tuple of nothing
        pos = new ArrayList<>();
        RCs = new ArrayList<>();
    }
    
    
    public RCTuple(ArrayList<Integer> pos, ArrayList<Integer> RCs) {
        this.pos = pos;
        this.RCs = RCs;
    }
    
    
    //one-RC tuple
    public RCTuple(int pos1, int RC1){
    	this();
    	set(pos1, RC1);
    }
    
    
    //a pair
    public RCTuple(int pos1, int RC1, int pos2, int RC2){
    	this();
    	set(pos1, RC1, pos2, RC2);
    }

    // a triple
	public RCTuple(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
    	this();
    	set(pos1, rc1, pos2, rc2, pos3, rc3);
	}

	// a quad
	public RCTuple(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, int pos4, int rc4) {
		this();
		set(pos1, rc1, pos2, rc2, pos3, rc3, pos4, rc4);
	}

    //Sometimes we'll want to generate an RC tuple from a conformation, specified as RCs for all positions
    //in order.  
    //In this case, negative values are not (fully) defined, so the tuple contains all positions
    //with positive values in conf
    public RCTuple(int[] conf){
    	this();
    	set(conf);
    }

    public RCTuple(ConfIndex index) {
    	this();
    	set(index);
	}
    
    public RCTuple set(int pos, int rc) {
    	this.pos.clear();
    	this.RCs.clear();
    	
        this.pos.add(pos);
        this.RCs.add(rc);

        return this;
    }
    
    public RCTuple set(int pos1, int rc1, int pos2, int rc2) {
    	this.pos.clear();
    	this.RCs.clear();
    	
        this.pos.add(pos1);
        this.RCs.add(rc1);
        
        this.pos.add(pos2);
        this.RCs.add(rc2);

        return this;
    }

	public RCTuple set(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
		this.pos.clear();
		this.RCs.clear();

		this.pos.add(pos1);
		this.RCs.add(rc1);

		this.pos.add(pos2);
		this.RCs.add(rc2);

		this.pos.add(pos3);
		this.RCs.add(rc3);

		return this;
	}

	public RCTuple set(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, int pos4, int rc4) {
		this.pos.clear();
		this.RCs.clear();

		this.pos.add(pos1);
		this.RCs.add(rc1);

		this.pos.add(pos2);
		this.RCs.add(rc2);

		this.pos.add(pos3);
		this.RCs.add(rc3);

		this.pos.add(pos4);
		this.RCs.add(rc4);

		return this;
	}

	public void set(int[] conf) {
    	pos.clear();
    	RCs.clear();
        for(int posNum=0; posNum<conf.length; posNum++){
            if(conf[posNum]>=0){//RC fully defined
                pos.add(posNum);
                RCs.add(conf[posNum]);
            }
        }
    }
    
    public void set(RCTuple other) {
    	pos.clear();
    	RCs.clear();
    	pos.addAll(other.pos);
    	RCs.addAll(other.RCs);
    }

    public void set(ConfIndex index) {
    	pos.clear();
    	RCs.clear();
    	for (int i=0; i<index.numDefined; i++) {
    		pos.add(index.definedPos[i]);
    		RCs.add(index.definedRCs[i]);
		}
	}
    
    public int size() {
    	return pos.size();
    }


	/**
	 * returns true if they are the same tuple AND if positions are in the same order
	 */
	public boolean isSameTuple(RCTuple tuple2){
    	
    	// short circuit: same instance must have same value
    	if (this == tuple2) {
    		return true;
    	}
    	
        //do the two tuple objects specify the same tuple of RCs?
        if( (pos.size()!=RCs.size()) || (tuple2.pos.size()!=tuple2.RCs.size()) )
            throw new RuntimeException("ERROR: Ill-defined RC tuple");
        
        
        if(pos.size()!=tuple2.pos.size())
            return false;
        
        //tuples are well-defined and same size...check position by position
        for(int index=0; index<pos.size(); index++){
            if(pos.get(index)!=tuple2.pos.get(index))
                return false;
            if(RCs.get(index)!=tuple2.RCs.get(index))
                return false;
        }
        
        //if we get here they're the same
        return true;
    }
    
    
    public String stringListing(){
        //Listing the RCs in a string
        String ans = "";
        
        for(int posNum=0; posNum<pos.size(); posNum++){
            ans = ans + "Res " + pos.get(posNum) + " RC " + RCs.get(posNum) + " ";
        }
        
        return ans;
    }
    
    
    public RCTuple subtractMember(int index){
        //Make a copy of this RCTuple with the given member removed
        //index is an index in pos and RCs
        ArrayList<Integer> newPos = new ArrayList<>();
        ArrayList<Integer> newRCs = new ArrayList<>();
        
        for(int ind=0; ind<pos.size(); ind++){
            if(ind!=index){
                newPos.add(pos.get(ind));
                newRCs.add(RCs.get(ind));
            }
        }
        
        return new RCTuple(newPos,newRCs);
    }
    
    @SuppressWarnings("unchecked")
	public RCTuple addRC(int addedPos, int addedRC){
        //Make a copy of this RCTuple with (addPos,addRC) added
        ArrayList<Integer> newPos = (ArrayList<Integer>) pos.clone();
        ArrayList<Integer> newRCs = (ArrayList<Integer>) RCs.clone();
        
        newPos.add(addedPos);
        newRCs.add(addedRC);
        
        return new RCTuple(newPos,newRCs);
    }

    @Override
	public int hashCode() {
		return HashCalculator.combineHashes(
			pos.hashCode(),
			RCs.hashCode()
		);
	}

	@Override
	public boolean equals(Object other) {
    	return other instanceof RCTuple && equals((RCTuple)other);
	}

	public boolean equals(RCTuple other) {
    	return isSameTuple(other);
	}

	@Override
	public String toString() {
		StringBuilder buf = new StringBuilder();
		buf.append("[");
		for (int i=0; i<size(); i++) {
			if (i > 0) {
				buf.append(",");
			}
			buf.append(pos.get(i));
			buf.append("=");
			buf.append(RCs.get(i));
		}
		buf.append("]");
		return buf.toString();
	}

	public void sortPositions() {

		// sort the positions using a simple insertion sort
		// tuples are always small (n << 100), so insertion sort should be fast enough
		// NOTE: we need to sort two arrays simultaneously, so we can't use any library sorts
		int n = size();
		for (int i=1; i<n; i++) {

			int tempPos = pos.get(i);
			int tempRC = RCs.get(i);

			int j;
			for (j=i; j>=1 && tempPos < pos.get(j-1); j--) {
				pos.set(j, pos.get(j-1));
				RCs.set(j, RCs.get(j-1));
			}
			pos.set(j, tempPos);
			RCs.set(j, tempRC);
		}
	}

	public RCTuple sorted() {
		sortPositions();
		return this;
	}

	public void checkSortedPositions() {
		for (int i=1; i<pos.size(); i++) {
			if (pos.get(i) <= pos.get(i - 1)) {
				throw new IllegalStateException("RCTuple positions are not sorted");
			}
		}
	}

	public RCTuple intersect(RCTuple other) {
		return RCTuple.intersect(this, other);
	}

	public static RCTuple intersect(RCTuple first, RCTuple second) {
		RCTuple out = new RCTuple();
		for(int tupIndex = 0; tupIndex < first.size(); tupIndex++)
		{
			int firstPos = first.pos.get(tupIndex);
			int firstRC  = first.RCs.get(tupIndex);
			if(second.pos.contains(firstPos) && second.RCs.get(tupIndex) == firstRC)
				out = out.addRC(firstPos, firstRC);
		}
		return out;
	}

	public Integer getRC(int index) {
		for (int i=0; i<size(); i++) {
			if (pos.get(i) == index) {
				return RCs.get(i);
			}
		}
		return null;
	}
}
