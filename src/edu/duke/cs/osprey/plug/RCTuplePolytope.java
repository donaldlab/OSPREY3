/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import static edu.duke.cs.osprey.plug.LPChecks.constrAsFunction;
import java.io.Serializable;
import java.util.ArrayList;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * Polytope for a particular set of DOFs, representing the feasible region for an RC tuple
 * 
 * @author mhall44
 */
public class RCTuplePolytope implements Serializable {
    
    ArrayList<DegreeOfFreedom> DOFs;
    ArrayList<LinearConstraint> constr;
    
    ArrayList<String> atomPairNames;//DEBUG!!!  names of atom pairs

    public RCTuplePolytope(ArrayList<DegreeOfFreedom> DOFs, ArrayList<LinearConstraint> constr) {
        this.DOFs = DOFs;
        this.constr = constr;
    }
    
    /*ArrayList<LinearConstraint> expandConstraints(ArrayList<DegreeOfFreedom> allDOFs){
        //express as constraints in a larger space of DOFs. 
        ArrayList<LinearConstraint> ans = new ArrayList<>();
        for(LinearConstraint c : constr){
            double bigGrad[] = new double[allDOFs.size()];
            for(int dof=0; dof<DOFs.size(); dof++)
                bigGrad[allDOFs.indexOf(DOFs.get(dof))] = c.getCoefficients().getEntry(dof);
            ans.add(new LinearConstraint(bigGrad, c.getRelationship(), c.getValue()));
        }
        return ans;
    }*/
    
    ArrayList<LinearConstraint> expandConstraints(ArrayList<DegreeOfFreedom> allDOFs){
        return transferConstraints(DOFs,allDOFs,constr);
    }
    
    public static ArrayList<LinearConstraint> transferConstraints(ArrayList<DegreeOfFreedom> oldDOFs,
            ArrayList<? extends DegreeOfFreedom> newDOFs, ArrayList<LinearConstraint> oldConstr){
        //transfer constraints from oldDOFs to newDOFs, applying 0's if new dof not among old dofs
        ArrayList<LinearConstraint> ans = new ArrayList<>();
        for(LinearConstraint c : oldConstr){
            double bigGrad[] = new double[newDOFs.size()];
            for(int dof=0; dof<oldDOFs.size(); dof++){
                /*int index = newDOFs.indexOf(oldDOFs.get(dof));
                if(index!=-1)
                    bigGrad[index] = c.getCoefficients().getEntry(dof);*/
                for(int index=0; index<newDOFs.size(); index++){//DEBUG!! compare by confDOFNum so can compare across confSpaces
                    if(oldDOFs.get(dof).getName().equalsIgnoreCase(newDOFs.get(index).getName())){
                        bigGrad[index] = c.getCoefficients().getEntry(dof);
                        break;
                    }
                }
            }
            ans.add(new LinearConstraint(bigGrad, c.getRelationship(), c.getValue()));
        }
        return ans;
    }
    
    
    public void checkDOFList(ArrayList<DegreeOfFreedom> checkDOFs){
        //check that DOFs matches check DOFs
        if(checkDOFs.size() != DOFs.size())
            throw new RuntimeException("ERROR: Wrong number of DOFs for RCTuplePolytope");
        for(int dof=0; dof<DOFs.size(); dof++){
            if(DOFs.get(dof)!=checkDOFs.get(dof))
                throw new RuntimeException("ERROR: DOFs don't match");
        }
    }
    
    
    public boolean containsPoint(double[] pt){
        //Does the polytope contain the point (specified as DOF values)?
        return LPChecks.isPointInPolytope(constr,pt);
    }
    
    
    public void listClashes(double[] pt){
        for(int c=0; c<constr.size(); c++){
            if(constrAsFunction(constr.get(c)).value(pt) > 1e-10)//inequality violated
                System.out.println("INEQUALITY VIOLATED: "+atomPairNames.get(c));
        }
    }

    public ArrayList<LinearConstraint> getConstr() {
        return constr;
    }
    
    public boolean hasNonBoxConstr(){
        //DEBUG!!  Assumes storing voxel constr here
        return constr.size() > 2*DOFs.size();
    }

    public ArrayList<String> getAtomPairNames() {
        return atomPairNames;
    }

    public ArrayList<DegreeOfFreedom> getDOFs() {
        return DOFs;
    }
    
    
    
}
