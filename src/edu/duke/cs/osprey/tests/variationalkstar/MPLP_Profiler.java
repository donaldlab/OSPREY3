/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests.variationalkstar;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.partitionfunctionbounds.GumbelMapTree;
import edu.duke.cs.osprey.partitionfunctionbounds.MapPerturbation;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author hmn5
 */
public class MPLP_Profiler {

    static EnergyMatrix emat;
    static PruningMatrix pruneMat;
    static ArrayList<ArrayList<Integer>> unprunedRCsAtPos = new ArrayList<>();
    static int numPos;

    public static void main(String[] args)
            throws Exception {

        String path = new File("").getAbsolutePath();
        if (!path.endsWith("MPLP/Profiling")) {
            throw new Error("This profiler was designed to be run in test/4HEM/Profiler folder\n\tcwd: " + path);
        }

        //load configurations
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();

        SearchProblem searchProb = cfp.getSearchProblem();

        double pruningInterval = 50;
        precomputeMatrices(searchProb, cfp, pruningInterval);
        emat = searchProb.emat;
        pruneMat = searchProb.pruneMat;
        numPos = emat.numPos();
        for (int pos = 0; pos < numPos; pos++) {
            unprunedRCsAtPos.add(pruneMat.unprunedRCsAtPos(pos));
        }

        ConfTree.mplpScore = true;
        ConfTree.traditionalScore = false;
        ConfTree tree = new ConfTree(searchProb);
        if (true) {
            Stopwatch.start();
            AStarNode rootNode = tree.rootNode();
            int[] blankConf = new int[searchProb.emat.numPos()];
            Arrays.fill(blankConf, -1);
            tree.initQueue(rootNode);
            Stopwatch.stop();
            System.out.println("MPLP Finished in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));

            Stopwatch.start();
            double tradScore = scoreConf(blankConf);
            Stopwatch.stop();
            System.out.println("Traditional Score: " + tradScore);
            System.out.println("Traditional Finished in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));
        }
        else if (false){
            Stopwatch.start();
            tree.nextConf();
            Stopwatch.stop();
            System.out.println("GMEC Found in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));
        } else {
            System.out.println("Starting Gumbel MAP Tree");
//            GumbelMapTree gTree = new GumbelMapTree(searchProb);
//            gTree.nextConf();

            cfp.params.setValue("STERICTHRESH", "1000");
            MapPerturbation mp = new MapPerturbation(searchProb);
            double logZ = mp.calcUBLogZ(20);
            System.out.println("UpperBound: "+logZ);
            
            PartFuncTree.verbose = true;
            PartFuncTree pfTree = new PartFuncTree(searchProb);
            double logZPF  = pfTree.computeEpsilonApprox(0.8);
            System.out.println("LogZ PartFunc: "+logZPF);
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

    static double scoreConf(int[] partialConf) {
        RCTuple definedTuple = new RCTuple(partialConf);

        double score = emat.getConstTerm() + emat.getInternalEnergy(definedTuple);//"g-score"

        //score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
        //plus contributions associated with each of the undefined res ("h-score")
        for (int level = 0; level < numPos; level++) {
            if (partialConf[level] < 0) {//level not fully defined

                double resContribLB = Double.POSITIVE_INFINITY;//lower bound on contribution of this residue
                //resContribLB will be the minimum_{rc} of the lower bound assuming rc assigned to this level

                for (int rc : unprunedRCsAtPos.get(level)) {
                    resContribLB = Math.min(resContribLB, RCContributionLB(level, rc, definedTuple, partialConf));
                }

                score += resContribLB;
            }
        }

        return score;
    }

    static double RCContributionLB(int level, int rc, RCTuple definedTuple, int[] partialConf) {
        //Provide a lower bound on what the given rc at the given level can contribute to the energy
        //assume partialConf and definedTuple

        double rcContrib = emat.getOneBody(level, rc);

        //for this kind of lower bound, we need to split up the energy into the defined-tuple energy
        //plus "contributions" for each undefined residue
        //so we'll say the "contribution" consists of any interactions that include that residue
        //but do not include higher-numbered undefined residues
        for (int level2 = 0; level2 < numPos; level2++) {

            if (partialConf[level2] >= 0 || level2 < level) {//lower-numbered or defined residues

                double levelBestE = Double.POSITIVE_INFINITY;//best pairwise energy
                ArrayList<Integer> allowedRCs = allowedRCsAtLevel(level2, partialConf);

                for (int rc2 : allowedRCs) {

                    double interactionE = emat.getPairwise(level, rc, level2, rc2);

                    double higherLB = higherOrderContribLB(partialConf, level, rc, level2, rc2);
                    //add higher-order terms that involve rc, rc2, and parts of partialConf

                    interactionE += higherLB;

                    //besides that only residues in definedTuple or levels below level2
                    levelBestE = Math.min(levelBestE, interactionE);
                }

                rcContrib += levelBestE;
            }
        }

        return rcContrib;
    }

    static double higherOrderContribLB(int[] partialConf, int pos1, int rc1, int pos2, int rc2) {
        //higher-order contribution for a given RC pair, when scoring a partial conf

        HigherTupleFinder<Double> htf = emat.getHigherOrderTerms(pos1, rc1, pos2, rc2);

        if (htf == null) {
            return 0;//no higher-order interactions
        } else {
            return higherOrderContribLB(partialConf, htf, pos2);
        }
    }

    static double higherOrderContribLB(int[] partialConf, HigherTupleFinder<Double> htf, int level2) {
        //recursive function to get lower bound on higher-than-pairwise terms
        //this is the contribution to the lower bound due to higher-order interactions
        //of the RC tuple corresponding to htf with "lower-numbered" residues (numbering as in scoreConf:
        //these are residues that are fully defined in partialConf, or are actually numbered <level2)

        double contrib = 0;

        for (int iPos : htf.getInteractingPos()) {//position has higher-order interaction with tup
            if (posComesBefore(iPos, level2, partialConf)) {//interaction in right order
                //(want to avoid double-counting)

                double levelBestE = Double.POSITIVE_INFINITY;//best value of contribution
                //from tup-iPos interaction
                ArrayList<Integer> allowedRCs = allowedRCsAtLevel(iPos, partialConf);

                for (int rc : allowedRCs) {

                    double interactionE = htf.getInteraction(iPos, rc);

                    //see if need to go up to highers order again...
                    HigherTupleFinder htf2 = htf.getHigherInteractions(iPos, rc);
                    if (htf2 != null) {
                        interactionE += higherOrderContribLB(partialConf, htf2, iPos);
                    }

                    //besides that only residues in definedTuple or levels below level2
                    levelBestE = Math.min(levelBestE, interactionE);
                }

                contrib += levelBestE;//add up contributions from different interacting positions iPos
            }
        }

        return contrib;
    }

    static ArrayList<Integer> allowedRCsAtLevel(int level, int[] partialConf) {
        //What RCs are allowed at the specified level (i.e., position num) in the given partial conf?
        ArrayList<Integer> allowedRCs;

        if (partialConf[level] == -1)//position undefined: consider all RCs
        {
            allowedRCs = unprunedRCsAtPos.get(level);
        } else if (partialConf[level] >= 0) {
            allowedRCs = new ArrayList<>();
            allowedRCs.add(partialConf[level]);
        } else {
            throw new UnsupportedOperationException("ERROR: Partially assigned position not yet supported in A*");
        }

        return allowedRCs;
    }

    static private boolean posComesBefore(int pos1, int pos2, int partialConf[]) {
        //for purposes of contributions to traditional conf score, 
        //we go through defined and then through undefined positions (in partialConf);
        //within each of these groups we go in order of position number
        if (partialConf[pos2] >= 0) {//pos2 defined
            return (pos1 < pos2 && partialConf[pos1] >= 0);//pos1 must be defined to come before pos2
        } else//pos1 comes before pos2 if it's defined, or if pos1<pos2
        {
            return (pos1 < pos2 || partialConf[pos1] >= 0);
        }
    }

}
