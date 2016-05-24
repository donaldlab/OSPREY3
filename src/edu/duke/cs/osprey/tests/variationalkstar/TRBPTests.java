/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests.variationalkstar;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.partitionfunctionbounds.MarkovRandomField;
import edu.duke.cs.osprey.partitionfunctionbounds.TRBP;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.Stopwatch;
import java.io.File;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author hmn5
 */
public class TRBPTests {

    public static void main(String[] args)
            throws Exception {

        String path = new File("").getAbsolutePath();
        if (!path.endsWith("VariationalKStar/TRBP")){
            throw new Error("This test was designed to be run in test/VariationalKStar/TRBP folder\n\tcwd: " + path);
        }

        //load configurations
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();

        SearchProblem searchProb = cfp.getSearchProblem();

        double pruningInterval = 50;
        precomputeMatrices(searchProb, cfp, pruningInterval);

        MarkovRandomField mrf = new MarkovRandomField(searchProb, 0.0);

        TRBP.setNumEdgeProbUpdates(0);
        TRBP.verbose = true;
        Stopwatch.start();
        TRBP trbp = new TRBP(mrf);
        Stopwatch.stop();

        double logZ = trbp.getLogZ();
        double correctLogZ = 906.6176159038234;

        System.out.println("Finished in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));
        System.out.println("LogZ: " + logZ);
        System.out.println("Error: " + Math.abs(correctLogZ - logZ));
        if (Math.abs(correctLogZ - logZ) < 0.1){
            System.out.println("Test Passed: Small Errors Results from different parameter settings");
        }
        else{
            System.out.println("Test Failed: This Could be an error or it could be a large change in parameters");
        }
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
