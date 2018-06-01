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

package edu.duke.cs.osprey.pruning;

/**
 *
 * @author mhall44
 */
public class PruningMethod {
    
    boolean useCompetitor;
    int numPos;//Are we pruning singles, pairs, triples, or what?
    CheckSumType cst;
    
    enum CheckSumType {
        GOLDSTEIN, INDIRECT, CONFSPLIT1, CONFSPLIT2, BOUNDS;    
    }

    public PruningMethod(boolean useCompetitor, int numPos, CheckSumType cst) {
        this.useCompetitor = useCompetitor;
        this.numPos = numPos;
        this.cst = cst;
    }
    
    
    
    
    public static PruningMethod getMethod(String name){
        //make PruningMethod object with settings implied by the name: "GoldsteinPairs", "Bounds", etc.
        
        if(name.equalsIgnoreCase("GOLDSTEIN"))
            return new PruningMethod(true,1,CheckSumType.GOLDSTEIN);
        
        else if(name.equalsIgnoreCase("GOLDSTEIN PAIRS"))
            return new PruningMethod(true,2,CheckSumType.GOLDSTEIN);
        
        else if(name.equalsIgnoreCase("GOLDSTEIN PAIRS FULL"))//SYNONYM FOR PAIRS
            return new PruningMethod(true,2,CheckSumType.GOLDSTEIN);
        
        else if(name.equalsIgnoreCase("GOLDSTEIN TRIPLES"))
            return new PruningMethod(true,3,CheckSumType.GOLDSTEIN);
        
        else if(name.equalsIgnoreCase("INDIRECT"))
            return new PruningMethod(true,1,CheckSumType.INDIRECT);
        
        else if(name.equalsIgnoreCase("INDIRECT PAIRS"))
            return new PruningMethod(true,2,CheckSumType.INDIRECT);
        
        else if(name.equalsIgnoreCase("SPLIT1"))
            return new PruningMethod(true,1,CheckSumType.CONFSPLIT1);
        
        else if(name.equalsIgnoreCase("SPLIT2"))
            return new PruningMethod(true,1,CheckSumType.CONFSPLIT2);
        
        else if(name.equalsIgnoreCase("BOUNDS"))
            return new PruningMethod(false,1,CheckSumType.BOUNDS);
        
        else if(name.equalsIgnoreCase("BOUNDING FLAGS"))//aka bounds pairs
            return new PruningMethod(false,2,CheckSumType.BOUNDS);
        
        throw new RuntimeException("ERROR: Unrecognized pruning method: "+name);
    }
    
    
    String name(){
        //a name for the pruning method
        return cst.name()+" FOR "+numPos+" POSITION(S)";
    }
    
}
