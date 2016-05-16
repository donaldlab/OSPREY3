/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.partitionfunctionbounds.MarkovRandomField;
import edu.duke.cs.osprey.partitionfunctionbounds.TRBP_Refactor_2;
import edu.duke.cs.osprey.partitionfunctionbounds.TRBP_Refactor_3;
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
        if (!path.endsWith("4HEM/TRBP_Test")) {
            throw new Error("This profiler was designed to be run in test/4HEM/TRBP_Test folder\n\tcwd: " + path);
        }

        //load configurations
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();

        SearchProblem searchProb = cfp.getSearchProblem();

        double pruningInterval = 50;
        precomputeMatrices(searchProb, cfp, pruningInterval);

        MarkovRandomField mrf = new MarkovRandomField(searchProb, 0.0);

        TRBP_Refactor_3.setNumEdgeProbUpdates(0);
        Stopwatch.start();
        TRBP_Refactor_3 trbp = new TRBP_Refactor_3(mrf);
        Stopwatch.stop();

        System.out.println("Finished in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));
        System.out.println("LogZ: "+trbp.getLogZ());
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
