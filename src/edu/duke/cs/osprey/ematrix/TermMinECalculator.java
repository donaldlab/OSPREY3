/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.handlempi.MPISlaveTask;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MolecEObjFunction;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

/**
 *
 * @author mhall44
 */

//This class is meant to precalculate the minimized energy of an energy term
//(e.g. intra+shell or pairwise)
//it should be suitable for calculating at a slave node

public class TermMinECalculator implements MPISlaveTask {
    
    ConfSpace confSpace;
    boolean doingEPIC;//doing EPIC fit instead of just minimum computation
    boolean doingIntra;//doing just intra energy (only relevant for one-body energies)
    //if false, one-body energies will be intra+shell
    
    int res[];//the residues to handle (can be either 1 or 2)
    
    EnergyFunction termE;//evaluates the energy for this term
    
    //So this class can be used to calculate 
    
    
    public TermMinECalculator(ConfSpace s, ArrayList<Residue> shellResidues, 
            boolean doEPIC, boolean doIntra, int... resToCalc){
        
        confSpace = s;
        doingEPIC = doEPIC;
        doingIntra = doIntra;
        res = resToCalc;
        
        Residue firstRes = confSpace.posFlex.get(res[0]).res;
        
        if(doingIntra)//intra energy only
            termE = EnvironmentVars.curEFcnGenerator.singleResEnergy(firstRes);
        else{
            
            if(res.length==1){//intra+shell
                termE =EnvironmentVars.curEFcnGenerator.intraAndShellEnergy(firstRes, shellResidues);
            }
            else if(res.length==2){//pairwise
                Residue secondRes = confSpace.posFlex.get(res[1]).res;
                termE = EnvironmentVars.curEFcnGenerator.resPairEnergy(firstRes, secondRes);
            }
            else{
                throw new RuntimeException("ERROR: Can only precompute energy for 1- and 2-body terms");
                //we are excluding intra-shell interactions throughout this version of OSPREY,
                //since they don't change and are time-consuming
            }
        }
    }
    

    @Override
    public Object doCalculation() {
        //This operation, which can be handled on a slave node once this TermMinECalculator is sent there,
        //calculates the specified energies
                
        if(!doingEPIC){
            if(res.length==1)
                return oneBodyMinEnergies();
            else if(res.length==2)
                return pairwiseMinEnergies();
        }
        /*
        else {
            if(res.length==1)
                return oneBodyPoly();
            else if(res.length==2)
                return pairwisePoly();
        }
        */
        throw new RuntimeException("ERROR: Calculation type not understood in TermMinECalculator");
    }
    
    
    public ArrayList<Double> oneBodyMinEnergies(){
        //list minimized one-body energies for all the RCs in res
        //the energy minimized can be just intra or intra+shell (decided in constructor)
        
        ArrayList<Double> ans = new ArrayList<>();
        ArrayList<RC> RCList = confSpace.posFlex.get(res[0]).RCs;
        
        for(int RCNum=0; RCNum<RCList.size(); RCNum++){
            RCTuple RCTup = new RCTuple(res[0],RCNum);
            double minE = calcMinE(RCTup);
            ans.add(minE);
        }
        
        return ans;
    }
    
    
    public ArrayList<ArrayList<Double>> pairwiseMinEnergies(){
        //list minimized one-body energies for all the RCs in res
        
        ArrayList<ArrayList<Double>> ans = new ArrayList<>();
        ArrayList<RC> RCList1 = confSpace.posFlex.get(res[0]).RCs;
        ArrayList<RC> RCList2 = confSpace.posFlex.get(res[1]).RCs;
        
        for(int firstRCNum=0; firstRCNum<RCList1.size(); firstRCNum++){
            
            ArrayList<Double> EList = new ArrayList<>();
            
            for(int secondRCNum=0; secondRCNum<RCList2.size(); secondRCNum++){
                RCTuple RCTup = new RCTuple(res[0],firstRCNum,res[1],secondRCNum);
                double minE = calcMinE(RCTup);
                EList.add(minE);
            }
            
            ans.add(EList);
        }
        
        return ans;
    }
    
    
    public double calcMinE(RCTuple RCs){
        //calculate minimum energy for an RC tuple

        
        MolecEObjFunction mof = new MolecEObjFunction(termE,confSpace,RCs);
        
        DoubleMatrix1D bestDOFVals;
        
        if(mof.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
            CCDMinimizer ccdMin = new CCDMinimizer(mof,true);
            bestDOFVals = ccdMin.minimize();
        }
        else//molecule is already in the right, rigid conformation
            bestDOFVals = DoubleFactory1D.dense.make(0);
        
        double energy = mof.getValue(bestDOFVals);
        
        return energy;
    }
    
}
