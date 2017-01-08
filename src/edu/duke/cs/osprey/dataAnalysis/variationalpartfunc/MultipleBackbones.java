/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dataAnalysis.variationalpartfunc;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.pruning.PruningControl;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.lang.ArrayUtils;

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
        unPruneAllRots(spBB);
        double[] dofVals = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
        for (int dofIndex = 0; dofIndex < dofVals.length; dofIndex++) {
            double[] dofValsToPrune = ArrayUtils.remove(dofVals, dofIndex);
            UpdatedPruningMatrix upm = pruneRotsWithPerturbation(spBB, dofValsToPrune);
            ConfTree tree = new ConfTree(spBB.emat, upm);
            int[] conf = tree.nextConf();
            if (conf != null) {
                System.out.println("Outputing GMEC For DOF " + dofVals[dofIndex]);
                spBB.outputMinimizedStruct(conf, spBB.name + "_" + dofVals[dofIndex] + ".GMEC.pdb");
            }
            else { 
                System.out.println("No GMEC Found For DOF "+ dofVals[dofIndex]);
            }
        }
    }

    //Prune all rotamers with a particular perturbation value;
    static UpdatedPruningMatrix pruneRotsWithPerturbation(SearchProblem sp, double[] pertVals) {
        UpdatedPruningMatrix upm = new UpdatedPruningMatrix(sp.pruneMat);
        for (int pos = 0; pos < sp.confSpace.numPos; pos++) {
            ArrayList<Integer> unprunedRCs = sp.pruneMat.unprunedRCsAtPos(pos);
            for (int rc : unprunedRCs) {
                for (double dof : sp.confSpace.posFlex.get(pos).RCs.get(rc).DOFmax) {
                    for (double pertVal : pertVals) {
                        if (dof == pertVal) {
                            upm.markAsPruned(new RCTuple(pos, rc));
                            break;
                        }
                    }
                }
            }
        }
        return upm;
    }

    static void unPruneAllRots(SearchProblem sp){
        for (int pos=0; pos<sp.confSpace.numPos; pos++){
            for (int rc=0; rc<sp.confSpace.posFlex.get(pos).RCs.size(); rc++){
                sp.pruneMat.setOneBody(pos, rc, Boolean.FALSE);
            }
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
