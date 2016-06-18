/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dataAnalysis.variationalpartfunc;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.io.File;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class MultipleBackbones {

    public static void main(String[] args)
            throws Exception {

        String path = new File("").getAbsolutePath();
/*        if (!path.endsWith("DEEPer/4HEM/Debug")) {
            throw new Error("This test was designed to be run in /usr/project/xtmp/hmn5/DEEPer/4HEM/Debug folder\n\tcwd: " + path);
        }*/

        ConfigFileParser cfpBB = setupRunWithDEEPer(path + "/");

        SearchProblem spBB = cfpBB.getSearchProblem();

        precomputeMatrices(spBB, cfpBB, Double.POSITIVE_INFINITY);
        Pruner deeBB = new Pruner(spBB, false, cfpBB.params.getDouble("BOUNDSTHRESH"),
                Double.POSITIVE_INFINITY, false, false);

        
        ConfTree tree = new ConfTree(spBB.emat, spBB.pruneMat);
        int[] conf;
        do {
            conf = tree.nextConf();
            if (conf != null) {
                double energy = spBB.emat.getInternalEnergy(new RCTuple(conf));
                System.out.println("Energy: " + energy);
            }
        } while (conf != null);

        return;
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
            PruningControl compPruning = cfp.setupPruning(searchSpace, 0, false, false);
            compPruning.setOnlyGoldstein(true);
            compPruning.prune();
            searchSpace.competitorPruneMat = searchSpace.pruneMat;
            searchSpace.pruneMat = null;
        }

        //Next, do DEE, which will fill in the pruning matrix
        PruningControl pruning = cfp.setupPruning(searchSpace, pruningInterval, false, false);

        pruning.prune();//pass in DEE options, and run the specified types of DEE            
    }

    private static ConfigFileParser setupRunWithDEEPer(String run) {
        String dee = run + "DEE.cfg";
        String sys = run + "System.cfg";
        String kstar = run + "KStar.cfg";

        String[] newArgs = {"-c", kstar, " ", dee, sys};
        ConfigFileParser cfp = new ConfigFileParser(newArgs);
        String pdbFile = cfp.params.getValue("PDBNAME");
        String pathToPDB = run + pdbFile;
        cfp.params.setValue("PDBNAME", pathToPDB);
        String name = cfp.params.getValue("RUNNAME");
        String pathName = run + name;
        cfp.params.setValue("RUNNAME", pathName);
        cfp.params.setValue("Ew", "1000000");
        cfp.loadData();
        return cfp;
    }

    private static ConfigFileParser setupRunWithoutDEEPer(String run) {
        String dee = run + "DEE.cfg";
        String sys = run + "System.cfg";
        String kstar = run + "KStar.cfg";

        String[] newArgs = {"-c", kstar, " ", dee, sys};
        ConfigFileParser cfp = new ConfigFileParser(newArgs);
        String pdbFile = cfp.params.getValue("PDBNAME");
        String pathToPDB = run + pdbFile;
        cfp.params.setValue("PDBNAME", pathToPDB);
        String name = cfp.params.getValue("RUNNAME");
        String pathName = run + name;
        cfp.params.setValue("RUNNAME", pathName);
        cfp.loadData();
        return cfp;
    }

}
