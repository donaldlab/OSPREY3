/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;
import java.util.ArrayList;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import java.io.Serializable;
import edu.duke.cs.osprey.tools.MinVolEllipse;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
/**
 *
 * @author hmn5
 */
public class SuperRC extends RC {
//This is a superRC class such that at each position we can have multiples residues
//The position will thus have a confSpace determined by SuperRCs 
//    
//I have implemented SuperRC and PositionConfSpaceSuper so that they each call 
//the normal PositionConfSpace and RC classes to create themselves.
//NOTE: Rather than extending RC it may be better to simply use composition or to have RCs subclass SuperRCs
    public ArrayList<String> AATypePerRes;
    public int[] rotNumPerRes;

    public ArrayList<ArrayList<DegreeOfFreedom>> DOFsPerRes;
    public ArrayList<ArrayList<Double>> DOFminPerRes;
    public ArrayList<ArrayList<Double>> DOFmaxPerRes;
    
    //Number of residues at position
    int numResiduesAtPosition;
    ArrayList<RC> rcList;
    
    public SuperRC(ArrayList<RC> rcList, int RCIndex){
        
        //RCIndex is the only field from RC that SuperRC relies on. Thus, if
        //we decide to separate the two, we just change RCIndex to SuperRCIndex 
        //and create a SuperRCIndex field
        this.RCIndex = RCIndex;
        this.rcList = rcList;
        this.numResiduesAtPosition = rcList.size();
        
        //Initialize the fields
        this.AATypePerRes = new ArrayList<String>(numResiduesAtPosition);
        this.rotNumPerRes = new int[numResiduesAtPosition];
        this.DOFsPerRes = new ArrayList<ArrayList<DegreeOfFreedom>>(numResiduesAtPosition);
        this.DOFminPerRes = new ArrayList<ArrayList<Double>>(numResiduesAtPosition);
        this.DOFmaxPerRes = new ArrayList<ArrayList<Double>>(numResiduesAtPosition);
        //Add info based on rcList
        for (int rcIndex = 0; rcIndex<rcList.size(); rcIndex++){
            RC rc = rcList.get(rcIndex);
            this.AATypePerRes.add(rc.AAType);
            this.rotNumPerRes[rcIndex] = rc.RCIndex;
            this.DOFsPerRes.add(rc.DOFs);
            this.DOFminPerRes.add(rc.DOFmin);
            this.DOFmaxPerRes.add(rc.DOFmax);
        }
    }
    
    public SuperRC(SuperRC rc1, SuperRC rc2, int RCIndex){
        ArrayList<String> AATypePerRes = new ArrayList<>();
        AATypePerRes.addAll(rc1.AATypePerRes);
        AATypePerRes.addAll(rc2.AATypePerRes);
        this.AATypePerRes = AATypePerRes;

        int[] rotNumPerRes = new int[this.AATypePerRes.size()];
        for(int rotNumIndex = 0; rotNumIndex<rc1.rotNumPerRes.length; rotNumIndex++){
           rotNumPerRes[rotNumIndex] = rc1.rotNumPerRes[rotNumIndex];
        }
        for(int rot2NumPerRes = 0; rot2NumPerRes<rc2.rotNumPerRes.length; rot2NumPerRes++){
            int rotNumIndex = rot2NumPerRes + rc1.rotNumPerRes.length;
            rotNumPerRes[rotNumIndex] = rc2.rotNumPerRes[rot2NumPerRes];
        }
        this.rotNumPerRes = rotNumPerRes;
        
        ArrayList<ArrayList<DegreeOfFreedom>> DOFsPerRes = new ArrayList<>();
        DOFsPerRes.addAll(rc1.DOFsPerRes);
        DOFsPerRes.addAll(rc2.DOFsPerRes);
        this.DOFsPerRes = DOFsPerRes;
        
        ArrayList<ArrayList<Double>> DOFminPerRes = new ArrayList<>();
        DOFminPerRes.addAll(rc1.DOFminPerRes);
        DOFminPerRes.addAll(rc2.DOFminPerRes);
        this.DOFminPerRes = DOFminPerRes;
        
        ArrayList<ArrayList<Double>> DOFmaxPerRes = new ArrayList<>();
        DOFmaxPerRes.addAll(rc1.DOFmaxPerRes);
        DOFmaxPerRes.addAll(rc2.DOFmaxPerRes);
        this.DOFmaxPerRes = DOFmaxPerRes;
        
        this.numResiduesAtPosition = rc1.numResiduesAtPosition + rc2.numResiduesAtPosition;
        
        ArrayList<RC> rcList = new ArrayList<>();
        rcList.addAll(rc1.rcList);
        rcList.addAll(rc2.rcList);
        this.rcList = rcList;
        
        this.RCIndex = RCIndex;
    }
    
    public boolean isParametricallyIncompatibleWith(SuperRC superRC2){
        //Two superRCs are parametrically incompatible if and only if there is a DOF that they share
        //for which they have different intervals
        //Thus, a conformation has a well-defined voxel if and only if it contains no parametrically
        //incompatible pairs of superRCs
        
        final double tol = 1e-8;

        //index over each rc in super-RC
        for (int rc1=0; rc1<this.DOFsPerRes.size(); rc1++){
            for (int rc2=0; rc2<rc1; rc2++){
                
                ArrayList<DegreeOfFreedom> DOFrc1 = this.DOFsPerRes.get(rc1);
                ArrayList<Double> DOFmin1 = this.DOFminPerRes.get(rc1);
                ArrayList<Double> DOFmax1 = this.DOFmaxPerRes.get(rc1);
                
                ArrayList<DegreeOfFreedom> DOFrc2 = superRC2.DOFsPerRes.get(rc2);
                ArrayList<Double> DOFmin2 = superRC2.DOFminPerRes.get(rc2);
                ArrayList<Double> DOFmax2 = superRC2.DOFmaxPerRes.get(rc2);
                
                for (int dof1=0; dof1<DOFrc1.size(); dof1++){
                    for (int dof2=0; dof2<DOFrc2.size(); dof2++){
                        if(DOFrc1.get(dof1) == DOFrc2.get(dof2)){// same DOF
                            
                            if ( Math.abs(DOFmin1.get(dof1) - DOFmin2.get(dof2)) > tol){
                                    return true;
                            }
                            if ( Math.abs(DOFmax1.get(dof1) - DOFmax2.get(dof2)) > tol){
                                return true;
                            }
                        }
                    }
                }
            }
        }
        //no incompatibility found!
        return false;
    }

}
