/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
