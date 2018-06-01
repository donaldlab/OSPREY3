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

package edu.duke.cs.osprey.plug;

import java.util.HashMap;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * This for checking for repeats of voxel constraints in a polytope
 * 
 * @author mhall44
 */
public class VoxConstrCache {
    
    //mapping DOF indices to their upper- and lower-bound values
    HashMap<Integer,Double> lb = new HashMap<>();
    HashMap<Integer,Double> ub = new HashMap<>();
    
    public VoxConstrCache(){}
    
    
    public boolean checkRedundancy(LinearConstraint constr){
        //check if the constraint is redundant with those already added
        //if not, and it is a voxel constraint, add it
        int voxConstrDOFNum = getVoxConstrDOFNum(constr.getCoefficients());
        if(voxConstrDOFNum==-1)//not a voxel constraint
            return false;
        switch(constr.getRelationship()){
            case GEQ:
                return checkRedundancy(voxConstrDOFNum, lb, constr.getValue());
            case LEQ:
                return checkRedundancy(voxConstrDOFNum, ub, constr.getValue());
            default:
                throw new RuntimeException("ERROR: Unexpected relationship");
        }
    }
    
    private static boolean checkRedundancy(int dofNum, HashMap<Integer,Double> curConstr, double val){
        //see if dofNum is listed in curConstr, and also look for contradictions
        if(curConstr.containsKey(dofNum)){
            if(Math.abs(val-curConstr.get(dofNum))>1e-10)
                throw new RuntimeException("ERROR: Conflicting voxel constraints");
            return true;
        }
        else{
            curConstr.put(dofNum, val);
            return false;
        }
    }
    
    private static int getVoxConstrDOFNum(RealVector coeff){
        //return what dof is constrained if this is a vox constr; else return -1
        int ans = -1;
        for(int dof=0; dof<coeff.getDimension(); dof++){
            double val = coeff.getEntry(dof);
            if(Math.abs(val)>1e-10){//nonzero
                if(ans!=-1)//this is not the first constrained dof.  Must not be a voxel constraint.  
                    return -1;
                if(Math.abs(val-1)>1e-10)//coeff 1 expected for vox constr
                    return -1;
                ans = dof;
            }
        }
        return ans;
    }
    
    
    
}
