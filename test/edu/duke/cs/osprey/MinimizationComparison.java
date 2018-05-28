/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;

/**
 *
 * Comparing minimization w/ and w/o EPIC
 * Set up to run for 1CC8.bbfree system
 * 
 * @author mhall44
 */
public class MinimizationComparison {
    
    
    public static void main(String args[]){
        //Trying to compute delta G's with continuous entropy between voxels.     
        //args like for findGMEC (currently set up for default 1CC8 system)
        
        MultiTermEnergyFunction.setNumThreads(4);
        
        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(args);//args are configuration files
        cfp.loadData();
        
        SearchProblem sp = cfp.getSearchProblem();
        sp.pruneMat = new PruningMatrix(sp.confSpace, Double.NEGATIVE_INFINITY);//no need to prune here
        sp.loadEnergyMatrix();
        sp.loadEPICMatrix();
        
        
        int conf1[] = new int[] {5,7,12,5,0,7,4};
        
        double Ereg = sp.minimizedEnergy(conf1);
        double Eepic = sp.EPICMinimizedEnergy(conf1);
        System.out.println("Regular E: "+Ereg);
        System.out.println("EPIC E: "+Eepic);
        
        
        //let's regenerate some of these energies, see if energy surfaces are actually different
        //or just minimized differently...
        RCTuple RCs = new RCTuple(conf1);
        MoleculeModifierAndScorer energy = new MoleculeModifierAndScorer(sp.fullConfE,sp.confSpace,RCs);
        Minimizer min = new CCDMinimizer(energy,false);
        DoubleMatrix1D optDOFVals = min.minimize().dofValues;
        double Ereg2 = energy.getValue(optDOFVals);
        System.out.println("Check regular E: "+Ereg2);
        
        EPICEnergyFunction epicEFunc = sp.epicMat.internalEnergyFunction(RCs, true);
        MoleculeModifierAndScorer epicEnergy = new MoleculeModifierAndScorer(epicEFunc,sp.epicMat.getConfSpace(),RCs);
        Minimizer epicMin = new CCDMinimizer(epicEnergy,false);
        DoubleMatrix1D epicOptDOFVals = epicMin.minimize().dofValues;
        double Eepic2 = epicEnergy.getValue(epicOptDOFVals);
        double EepicRegCheck = epicEnergy.getValue(optDOFVals);
        
        System.out.println("Check EPIC E: "+Eepic2);
        System.out.println("Check EPIC E at regular min: "+EepicRegCheck);
        
        ArrayList<Double> epicTermVals = epicEFunc.allTermValues();
        
        //Compare energies term by term...
        int numPosInTuple = RCs.pos.size();
        
        double Ereg3 = 0;
        double Eepic3 = 0;
        int termCount = 0;
        
        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            System.out.print("INTRA "+indexInTuple);
            
            int posNum = RCs.pos.get(indexInTuple);
            int RCNum = RCs.RCs.get(indexInTuple);
            
            //EPoly intraEPIC = sp.epicMat.getOneBody(posNum,RCNum);
            //double termEPIC = intraEPIC.evaluate(optDOFVals, true, true);
            double termEPIC = epicTermVals.get(termCount) + sp.emat.getOneBody(posNum, RCNum);
            termCount++;
            System.out.print(" EPIC E: "+termEPIC);
            Eepic3 += termEPIC;
            EnergyFunction intraE = EnvironmentVars.curEFcnGenerator.intraAndShellEnergy(sp.confSpace.posFlex.get(posNum).res, sp.shellResidues);
            double termReg = intraE.getEnergy();
            System.out.println(" ACTUAL E: "+termReg);
            Ereg3 += termReg;
            
            for(int index2=0; index2<indexInTuple; index2++){
                
                System.out.print("PAIRWISE "+indexInTuple+" " + index2);
                
                int pos2 = RCs.pos.get(index2);
                int rc2 = RCs.RCs.get(index2);
                
                //EPoly pairwiseEPIC = sp.epicMat.getPairwise(posNum,RCNum,pos2,rc2);
                //termEPIC = pairwiseEPIC.evaluate(optDOFVals, true, true);
                termEPIC = epicTermVals.get(termCount) + sp.emat.getPairwise(posNum, RCNum, pos2, rc2);
                termCount++;
                System.out.print(" EPIC E: "+termEPIC);
                Eepic3 += termEPIC;
                EnergyFunction pairwiseE = EnvironmentVars.curEFcnGenerator.resPairEnergy(sp.confSpace.posFlex.get(posNum).res, sp.confSpace.posFlex.get(pos2).res);
                termReg = pairwiseE.getEnergy();
                System.out.println(" ACTUAL E: "+termReg);
                Ereg3 += termReg;
            }
        }
        
        System.out.println("TOTAL ACTUAL E: "+Ereg3);
        System.out.println("TOTAL EPIC E: "+Eepic3);
    }
    
    
}
