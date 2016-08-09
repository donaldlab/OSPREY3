/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.GMECFinder;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.voxq.VoxelsDeltaG;

/**
 *
 * @author mhall44
 */
public class VoxDeltaGPlayground {
    
    
    public static void main(String args[]){
        //Trying to compute delta G's with continuous entropy between voxels.     
        //args like for findGMEC (currently set up for default 1CC8 system)
        
        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
	cfp.loadData();
        
        //Make two copies of the search problem so we can construct independent MoleculeModifierAndScorers
        //for each voxel
        SearchProblem sp1 = cfp.getSearchProblem();
        SearchProblem sp2 = cfp.getSearchProblem();
        
        System.out.println("Testing IVS...");
        int conf1[] = new int[] {5,7,12,5,1,7,4};
        int conf2[] = new int[] {5,7,7,5,0,7,4};
        //int conf2[] = new int[] {5,7,10,5,1,7,4};
        
        MoleculeModifierAndScorer mms1 = new MoleculeModifierAndScorer(sp1.fullConfE,
            sp1.confSpace, new RCTuple(conf1) );

        MoleculeModifierAndScorer mms2 = new MoleculeModifierAndScorer(sp2.fullConfE,
            sp2.confSpace, new RCTuple(conf2) );

        VoxelsDeltaG vdg = new VoxelsDeltaG(mms1,mms2);
        double dG = vdg.estDeltaG(0.05);
        System.out.println("delta G: "+dG);
        /*IntraVoxelSampler ivs = new IntraVoxelSampler(mms);
        for(int s=0; s<20; s++){
            System.out.println(ivs.nextSample());
        }*/
        System.exit(0);


        //conf2 = new int[] {5,7,12,5,1,7,4};
        //conf2 = new int[] {5,7,7,5,0,7,4};
    }
    
    
}
