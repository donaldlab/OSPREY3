/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.IdealSeparableReference;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.voxq.VoxelsDeltaG;
import java.util.Arrays;

/**
 *
 * Second attempt at BAR delta G's with LUTE. 
 * Uses regular (non-delta-G) pruning.  
 * 
 * @author mhall44
 */
public class ConfGTupleExpander2 extends ConfETupleExpander {
    
    EPICMatrix epicMat1, epicMat2;
    
        
    public ConfGTupleExpander2(SearchProblem sp){
        super(sp);
        
        epicMat1 = sp.epicMat;
        epicMat2 = (EPICMatrix) ObjectIO.deepCopy(epicMat1);
    }
    
    @Override
    double scoreAssignmentList(int[] assignmentList) {
        
        if(!sp.useEPIC)
            throw new RuntimeException("ERROR NEED EPIC FOR BAR+LUTE");
        
        MoleculeModifierAndScorer mms1 = new MoleculeModifierAndScorer(
                epicMat1.internalEnergyFunction(new RCTuple(assignmentList)), 
                epicMat1.getConfSpace(), new RCTuple(assignmentList) );
        
        CCDMinimizer ccdMin = new CCDMinimizer(mms1,false);
        DoubleMatrix1D center = ccdMin.minimize();
        MoleculeModifierAndScorer mms2 = new IdealSeparableReference(
                epicMat2.internalEnergyFunction(new RCTuple(assignmentList)), 
                epicMat2.getConfSpace(), new RCTuple(assignmentList), center );
        
        
        VoxelsDeltaG vdg = new VoxelsDeltaG(mms2,mms1,true);
        double E = vdg.estDeltaG(0.05);
        E += ((IdealSeparableReference)mms2).calcG();
        
        E += sp.emat.confE(assignmentList);//discrete part of energy

        if(E==Double.POSITIVE_INFINITY){//this is going to be a problem if used as a true value
            RCTuple tup = new RCTuple(assignmentList);
            if(isPruned(tup))
                throw new RuntimeException("ERROR: Scoring pruned conformation: "+tup.stringListing());
            else
                throw new RuntimeException("ERROR: Infinite E for unpruned conf: "+tup.stringListing());
        }

        return E;
    }


}

