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

package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import java.util.ArrayList;

/**
 * Implementation of a tuple expander that expands conformational energies--
 * including continuous flexibility, non-pairwise terms, etc.
 * 
 * @author mhall44
 */
public class ConfETupleExpander extends TupleExpander {
    
    SearchProblem sp;

    boolean fullPLUGPruning;//only allow sampled conformations that PLUG considers feasible


    public ConfETupleExpander(SearchProblem sp){
        super(sp.confSpace.numPos, sp.confSpace.getNumRCsAtPos(), sp.pruneMat.getPruningInterval(), sp.luteSettings);
        this.sp = sp;
        fullPLUGPruning = sp.plugMat!=null;//DEBUG!!!
        if(fullPLUGPruning)
            canCheckPartialPruning = false;
    }

    double worstELBDiff = 0;
    
    @Override
    double scoreAssignmentList(int[] assignmentList) {
        if(sp.useEPIC){//Faster if we can score by EPIC
            
            double E = sp.EPICMinimizedEnergy(assignmentList);

            RCTuple tup = new RCTuple(assignmentList);
            if(E==Double.POSITIVE_INFINITY){//this is going to be a problem if used as a true value
                if(isPruned(tup))
                    throw new RuntimeException("ERROR: Scoring pruned conformation: "+tup.stringListing());
                else
                    throw new RuntimeException("ERROR: Infinite E for unpruned conf: "+tup.stringListing());
            }

            if(fullPLUGPruning){
                if( ! sp.plugMat.isTupleFeasible(tup) ){
                    throw new RuntimeException("ERROR: Scoring PLUG-infeasible conformation: "+tup.stringListing());
                }
            }
            
            
            /*
            //code for debugging inaccurate pairwise lower bounds
            double LB = sp.lowerBound(assignmentList);
            if(E - LB  < -10){
                System.out.println("ERROR: E="+E+" LB="+LB);
                System.out.print( "Assignments: ");
                for(int a : assignmentList)
                    System.out.print(a+", ");
                
                RCTuple tup = new RCTuple(assignmentList);
                EPICEnergyFunction efunc = sp.epicMat.internalEnergyFunction(tup);
                MoleculeModifierAndScorer objFcn = new MoleculeModifierAndScorer(efunc,sp.epicMat.getConfSpace(),tup);
                
                Minimizer minim = new CCDMinimizer(objFcn,false);
                DoubleMatrix1D bestDOFVals = minim.minimize();
                System.out.println("Best DOF vals check: "+bestDOFVals);
                double Echeck = objFcn.getValue(bestDOFVals);
                System.out.println("E check: "+Echeck);
                
                System.out.println("EPIC term values: ");
                efunc.printAllTermValues();
                System.out.println("End EPIC term values");
                
                if(E-LB < worstELBDiff){
                    worstELBDiff = E-LB;
                    System.out.println("Outputting obj fcn to problem_obj_fcn"+worstELBDiff+".dat");
                    ObjectIO.writeObject(objFcn, "problem_obj_fcn"+worstELBDiff+".dat");
                }
            }*/
            
            return E;
        }
        else//Score by minimizing energy function directly
            return sp.minimizedEnergy(assignmentList);
    }

    
    @Override
    public boolean isPruned(RCTuple tup) {
        if (sp.pruneMat.isPruned(tup))
            return true;
        if(fullPLUGPruning){//pruning without marking--we only need expansion to be right on PLUG-feasible voxels
            if(!sp.plugMat.isTupleFeasible(tup))
                return true;
        }
        return false;
    }

    
    @Override
    void pruneTuple(RCTuple tup) {
        sp.pruneMat.markAsPruned(tup);
    }

    
    @Override
    ArrayList<RCTuple> higherOrderPrunedTuples(RCTuple tup) {
        //list higher-order pruned tuples containing the pair tup
        
        if(tup.pos.size() != 2)
            throw new RuntimeException("ERROR: higherOrderPrunedTuples is meant to take an RC pair as argument");
                
        HigherTupleFinder<Boolean> htf = 
                sp.pruneMat.getHigherOrderTerms(tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1));
        
        if(htf!=null){
            ArrayList<RCTuple> otherTups = htf.listInteractionsWithValue(true);
            //otherTups are recorded as what tup interacts with...add in tup to get the whole pruned tuple
            for(RCTuple otherTup: otherTups){
                for(int i=0; i<2; i++){
                    otherTup.pos.add(tup.pos.get(i));
                    otherTup.RCs.add(tup.RCs.get(i));
                }
            }
            
            return otherTups;
        }
        else//no interactions
            return new ArrayList<>();
    }
    
}
