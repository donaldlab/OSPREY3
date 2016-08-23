/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.voxq.VoxelsDeltaG;
import java.util.Arrays;

/**
 *
 * First attempt at BAR delta G's with LUTE. 
 * Uses regular (non-delta-G) pruning.  
 * 
 * @author mhall44
 */
public class ConfGTupleExpander1 extends ConfETupleExpander {
    
    EPICMatrix epicMat1, epicMat2;
    int[] refConf = null;//reference conformation
    
    
    boolean Gon = false;//init GMEC search is without G on
    
    public ConfGTupleExpander1(SearchProblem sp){
        super(sp);
        
        epicMat1 = sp.epicMat;
        epicMat2 = (EPICMatrix) ObjectIO.deepCopy(epicMat1);
    }
    
    @Override
    double scoreAssignmentList(int[] assignmentList) {
        
        if(!Gon)
            return super.scoreAssignmentList(assignmentList);
        
        if(!sp.useEPIC)
            throw new RuntimeException("ERROR NEED EPIC FOR BAR+LUTE");
            
        
        //OK FOR NOW ENERGY WILL BE JUST DELTA G RELATIVE TO REFERENCE CONF

                
        /*MoleculeModifierAndScorer mms1 = new MoleculeModifierAndScorer(sp1.fullConfE,
            sp1.confSpace, new RCTuple(conf1) );*/
        
        MoleculeModifierAndScorer mms1 = new MoleculeModifierAndScorer(
                epicMat1.internalEnergyFunction(new RCTuple(refConf)), 
                epicMat1.getConfSpace(), new RCTuple(refConf) );
        
        MoleculeModifierAndScorer mms2 = new MoleculeModifierAndScorer(
                epicMat2.internalEnergyFunction(new RCTuple(assignmentList)), 
                epicMat2.getConfSpace(), new RCTuple(assignmentList) );

        
        VoxelsDeltaG vdg = new VoxelsDeltaG(mms1,mms2,true);
        double E = vdg.estDeltaG(0.05);
        

        if(E==Double.POSITIVE_INFINITY){//this is going to be a problem if used as a true value
            RCTuple tup = new RCTuple(assignmentList);
            if(isPruned(tup))
                throw new RuntimeException("ERROR: Scoring pruned conformation: "+tup.stringListing());
            else
                throw new RuntimeException("ERROR: Infinite E for unpruned conf: "+tup.stringListing());
        }


        return E;
    }


    @Override
    double computeInitGMECEst(){//let's find this by some random iterations...
        double ans = Double.POSITIVE_INFINITY;
        TESampleSet tss = new TESampleSet(this);
        
        System.out.println("Computing initial GMEC estimate...");
        
        
        //Let's make sure, by DFS, that there is at least one unpruned conf
        //If not then our init GMEC is infinity
        int testSamp[] = new int[numPos];
        Arrays.fill(testSamp, -1);
        if( tss.finishSampleDFS(testSamp) == null )
            return Double.POSITIVE_INFINITY;
        
        
        for(int iter=0; iter<500/*0*/; iter++){
            int sample[] = new int[numPos];
            boolean success;
            do {
                Arrays.fill(sample,-1);
                success = tss.finishSample(sample);
            } while (!success);
            
            double score = scoreAssignmentList(sample);
            ans = Math.min(ans,score);
            
            
            //CHANGE IS HERE
            if(ans==score)//best so far
                refConf = sample;
        }
        
        //AND HERE
        Gon = true;
        
        return ans;
    }








}

