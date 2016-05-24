/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests.variationalkstar;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.partitionfunctionbounds.ReparamMRF;
import edu.duke.cs.osprey.partitionfunctionbounds.SCMF_Clamp;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.Stopwatch;
import java.io.File;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author hmn5
 */
public class SCMFTests {

    public static void main(String[] args)
            throws Exception {

        String path = new File("").getAbsolutePath();
        if (!path.endsWith("VariationalKStar/SCMF")){
            throw new Error("This tests was designed to be run in test/VariationalKStar/SCMF folder\n\tcwd: " + path);
        }

        //load configurations
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();

        SearchProblem searchProb = cfp.getSearchProblem();

        double pruningInterval = 50;
        precomputeMatrices(searchProb, cfp, pruningInterval);

        ReparamMRF mrf = new ReparamMRF(searchProb.emat, searchProb.pruneMat, 0.0);

        Stopwatch.start();
        SCMF_Clamp scmf = new SCMF_Clamp(mrf);
        Stopwatch.stop();

        double logZ = scmf.getLogZLB();
        double correctLogZ = 898.9578799814087;

        System.out.println("Finished Test 1 in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));
        System.out.println("LogZ: " + logZ);
        System.out.println("Error: " + Math.abs(correctLogZ - logZ));

        ReparamMRF mrf2 = new ReparamMRF(searchProb.emat, searchProb.pruneMat, 0.1);
        Stopwatch.start();
        SCMF_Clamp scmf2 = new SCMF_Clamp(mrf2);
        Stopwatch.stop();
        double logZ2 = scmf2.getLogZLB();
        double correctLogZ2 = 894.6048097520684;
        
        System.out.println("Finished Test 2 in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));
        System.out.println("LogZ: " + logZ2);
        System.out.println("Error: " + Math.abs(correctLogZ2 - logZ2));

    }

    static void precomputeMatrices(SearchProblem searchSpace, ConfigFileParser cfp, double pruningInterval) {
            //Precalculate TupleMatrices needed for GMEC computation.  Some of these may already be computed.  
        //All of these matrices except the basic pairwise energy matrix are pruning-dependent:
        //we can prune conformations whose energies are within pruningInterval
        //of the lowest pairwise lower bound

        //First calculate the pairwise energy matrix, if not already present
        searchSpace.loadEnergyMatrix();

        //Doing competitor pruning now
        //will limit us to a smaller, but effective, set of competitors in all future DEE
        if (searchSpace.competitorPruneMat == null) {
            System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
            PruningControl compPruning = cfp.setupPruning(searchSpace, 0, false, false);
            compPruning.setOnlyGoldstein(true);
            compPruning.prune();
            searchSpace.competitorPruneMat = searchSpace.pruneMat;
            searchSpace.pruneMat = null;
            System.out.println("COMPETITOR PRUNING DONE");
        }

        //Next, do DEE, which will fill in the pruning matrix
        PruningControl pruning = cfp.setupPruning(searchSpace, pruningInterval, false, false);

        pruning.prune();//pass in DEE options, and run the specified types of DEE            
    }
}
