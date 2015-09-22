/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.EllipseCoordDOF;
import edu.duke.cs.osprey.dof.MoveableStrand;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.CartesianProduct;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class PositionConfSpaceSuper extends PositionConfSpace{
    //Extension of PositionConfSpace to handle super-RCs and super-residues
    //Currently this class calls PositionConfSpace, but I think it could simply
    //replace PositionConfSpace
    public ArrayList<SuperRC> superRCs;
    public ArrayList<Residue> resList;//The list of residues involved
    
      
    //DOFIndices is indexed the same as resList but it returns an index into the 
    //confSpaceDOF
    public ArrayList<Integer> DOFIndices;
    ArrayList<ArrayList<EllipseCoordDOF>> ellipsoidalDOFsPerPos;
    boolean useEllipses = false;
    
    
    public PositionConfSpaceSuper(
    		ArrayList<Residue> resList,//perRes
    		ArrayList<ArrayList<DegreeOfFreedom>> resDOFsList,//perRes
    		ArrayList<ArrayList<String>> allowedAAsList,//perRes
                ArrayList<Integer> DOFIndices,//perRes
    		boolean contSCFlex,
                ArrayList<DegreeOfFreedom> strandDOFs,
                ArrayList<Perturbation> perts,
                ArrayList<ArrayList<double[]>> pertIntervals,
                ArrayList<ArrayList<ArrayList<int[]>>> pertStatesList,//perRes
                ArrayList<BBFreeBlock> bfbList,//perRes
    		boolean useEllipses) {
        
        
    	this.resList = resList;
        this.DOFIndices = DOFIndices;
        
        ArrayList<PositionConfSpace> positionConfSpacePerResidue = new ArrayList<>();

        ArrayList<ArrayList<RC>> rcSetPerRes = new ArrayList<ArrayList<RC>>();
        for (int resIndex = 0; resIndex < resList.size(); resIndex++){
            Residue res = resList.get(resIndex);
            ArrayList<DegreeOfFreedom> resDOFs = resDOFsList.get(resIndex);
            ArrayList<String> allowedAAs = allowedAAsList.get(resIndex);
            ArrayList<ArrayList<int[]>> pertStates = pertStatesList.get(resIndex);
            BBFreeBlock bfb = bfbList.get(resIndex);
            
            PositionConfSpace residueConfSpace = new PositionConfSpace(res, resDOFs, allowedAAs, contSCFlex, strandDOFs, perts, pertIntervals, pertStates, bfb, useEllipses);
            positionConfSpacePerResidue.add(residueConfSpace);
            if (useEllipses){
                this.useEllipses = true;
                ellipsoidalDOFsPerPos.add(residueConfSpace.ellipsoidalDOFs);
            }
            rcSetPerRes.add(positionConfSpacePerResidue.get(resIndex).RCs);
        }
        this.superRCs = new ArrayList<>();
        if (resList.size()==1){//If we just have one res per pos the each superRC is a list of one RC
            for(int superRCIndex=0;superRCIndex<rcSetPerRes.get(0).size(); superRCIndex++){
                ArrayList<RC> rcList = new ArrayList<>();
                rcList.add(rcSetPerRes.get(0).get(superRCIndex));
                SuperRC superRC = new SuperRC(rcList,superRCIndex);
                superRCs.add(superRC);
            }
        }
        else{//If not then we want the cartesian product of the RCs at each residue
            ArrayList<ArrayList<RC>> superRCList =  CartesianProduct.cartesianProduct(rcSetPerRes);
            for(int superRCIndex=0; superRCIndex<superRCList.size(); superRCIndex++){
                ArrayList<RC> rcList = superRCList.get(superRCIndex);
                SuperRC superRC = new SuperRC(rcList,superRCIndex);
                superRCs.add(superRC);
            }
        }
    }
    
    //This constructor is for merging position
    public PositionConfSpaceSuper(PositionConfSpaceSuper pcs1, PositionConfSpaceSuper pcs2){
        ArrayList<Residue> newResList = new ArrayList<>();
        newResList.addAll(pcs1.resList);
        newResList.addAll(pcs2.resList);
        this.resList = newResList;
        
        ArrayList<Integer> newDOFIndeces = new ArrayList<>();
        newDOFIndeces.addAll(pcs1.DOFIndices);
        newDOFIndeces.addAll(pcs2.DOFIndices);
        this.DOFIndices = newDOFIndeces;

        if (useEllipses){
            ArrayList<ArrayList<EllipseCoordDOF>> newEllipsoidalDOFsPerPos = new ArrayList<>();
            newEllipsoidalDOFsPerPos.addAll(pcs1.ellipsoidalDOFsPerPos);
            newEllipsoidalDOFsPerPos.addAll(pcs2.ellipsoidalDOFsPerPos);
        }
        
        //Now we create and add new super-RCs
        ArrayList<SuperRC> newSuperRCs = new ArrayList<>();
        int rcCounter = 0;
        for (int rcIndex1=0; rcIndex1<pcs1.superRCs.size(); rcIndex1++){
            SuperRC superRC1 = pcs1.superRCs.get(rcIndex1);
            for (int rcIndex2=0; rcIndex2<pcs2.superRCs.size(); rcIndex2++){
                SuperRC superRC2 = pcs2.superRCs.get(rcIndex2);
                if (!superRC1.isParametricallyIncompatibleWith(superRC2)) {
                    SuperRC newSuperRC = new SuperRC(superRC1, superRC2, rcCounter);
                    newSuperRCs.add(newSuperRC);
                    rcCounter++;
                }
            }
        }
        this.superRCs = newSuperRCs;
    }

    public ArrayList<ArrayList<EllipseCoordDOF>> getEllipsoidalArrayPerPos(){
        return this.ellipsoidalDOFsPerPos;
    }
}
