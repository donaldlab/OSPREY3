/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.util.ArrayList;

import edu.duke.cs.osprey.astar.ConfTreeSuper;
import edu.duke.cs.osprey.confspace.SearchProblemSuper;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.pruning.PruningControlSuper;
import edu.duke.cs.osprey.confspace.PositionConfSpaceSuper;
import edu.duke.cs.osprey.confspace.ConfSpaceSuper;
import edu.duke.cs.osprey.confspace.SuperRCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.partitionfunctionbounds.MarkovRandomField;
import edu.duke.cs.osprey.partitionfunctionbounds.SelfConsistentMeanField;
import edu.duke.cs.osprey.partitionfunctionbounds.GumbelDistribution;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.partitionfunctionbounds.MapPerturbation;
import edu.duke.cs.osprey.partitionfunctionbounds.SelfConsistentMeanField_Parallel;
import java.awt.Point;
import java.io.PrintStream;
import java.io.FileOutputStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.stream.Collectors;
import java.util.List;
import javafx.scene.input.KeyCode;

/**
 * KaDEEFinder computes the sequence with the highest K* score.using
 * probabilistic bounds.
 *
 * @author hmn5
 */
public class KaDEEFinder {

    ConfigFileParser cfp;

    // A search problem for Super RCs.
    SearchProblemSuper searchSpace;

    double Ew; // energy window for enumerating conformations: 0 for just GMEC
    double I0 = 0; // initial value of iMinDEE pruning interval
    boolean doIMinDEE;//do iMinDEE

    boolean useContFlex;
    //boolean enumByPairwiseLowerBound;//are we using a search method that 
    //enumerates by pairwise lower bound (minDEE-style)
    //or by (possibly approximated) true energy (rigid, EPIC, or tuple expander)?

    boolean outputGMECStruct;//write the GMEC structure to a PDB file

    boolean useEPIC = false;
    boolean useTupExp = false;

    boolean checkApproxE = true;//Use the actual energy function to evaluate
    //each enumerated conformation rather than just relying on the EPIC or tup-exp approximation

    boolean useEllipses = false;

    ExpFunction ef = new ExpFunction();
    double constRT = PoissonBoltzmannEnergy.constRT;

    public KaDEEFinder(ConfigFileParser cfp) {
        this.cfp = cfp;
        Ew = cfp.params.getDouble("Ew", 0);
        doIMinDEE = cfp.params.getBool("imindee", false);
        if (doIMinDEE) {
            I0 = cfp.params.getDouble("Ival", 5);
        }

        useContFlex = cfp.params.getBool("doMinimize", false);
        useTupExp = cfp.params.getBool("UseTupExp", false);
        useEPIC = cfp.params.getBool("UseEPIC", false);

        checkApproxE = cfp.params.getBool("CheckApproxE", true);

        if (doIMinDEE && !useContFlex) {
            throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility");
        }

        outputGMECStruct = cfp.params.getBool("OUTPUTGMECSTRUCT", false);

        useEllipses = cfp.params.getBool("useEllipses", false);
    }

    /**
     * Run KaDEE and compute the sequence with the highest K* score. TODO: For
     * now this class only computes the partition function for the sequence with
     * the highest K* score.
     */
    void doKaDEE() {
        double curInterval = I0;//For iMinDEE.  curInterval will need to be an upper bound

        searchSpace = cfp.getSearchProblemSuper();
        ConfSpaceSuper confSpaceSuper = searchSpace.confSpaceSuper;
        searchSpace.loadEnergyMatrix();
        for (int i = 0; i < 2; i++) {

            double pruningInterval = Double.POSITIVE_INFINITY;
            //Doing competitor pruning now
            //will limit us to a smaller, but effective, set of competitors in all future DEE
            if (searchSpace.competitorPruneMat == null) {
                System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
                PruningControlSuper compPruning = cfp.setupPruning(searchSpace, 0, false, false);
                compPruning.setOnlyGoldstein(true);
                compPruning.prune();
                searchSpace.competitorPruneMat = searchSpace.pruneMat;
                searchSpace.pruneMat = null;
                System.out.println("COMPETITOR PRUNING DONE");
            }
            //Next, do DEE, which will fill in the pruning matrix
            PruningControlSuper pruning = cfp.setupPruning(searchSpace, pruningInterval, false, false);
            pruning.prune();//pass in DEE options, and run the specified types of DEE 
//        BigInteger confSpace = new BigInteger("1");
//        for (int pos = 0; pos < searchSpace.emat.numPos(); pos++){
//            confSpace = confSpace.multiply(new BigInteger(((Integer) searchSpace.emat.oneBody.get(pos).size()).toString()));
//        }
            //System.out.println(confSpace);
            
            double gmecE = calcGMEC(searchSpace);
            System.out.println("GMEC Energy = "+gmecE);
            //Calculate GMEC (0.43429 = log_10(e) puts gmecE in log_10 scale 
            double gmecScore = -(0.4342944819) * gmecE/this.constRT;
            System.out.println("GMEC Bound = " + gmecScore);

            //SCMF
            double pruningInterval2 = 30.0;
            PruningControlSuper pruning2 = cfp.setupPruning(searchSpace, pruningInterval2, false, false);
            pruning2.prune();//pass in DEE options, and run the specified types of DEE 

            MarkovRandomField mrf = new MarkovRandomField(searchSpace, 0.0);
            SelfConsistentMeanField scmf = new SelfConsistentMeanField(mrf);
            scmf.run();
            BigDecimal Z = scmf.calcPartitionFunction();
            double logZLB = ef.log10(Z);

            System.out.println("Lower bound on log partition function (SCMF) = " + logZLB);
            //Next, do DEE, which will fill in the pruning matrix

            double pruningInterval3 = Double.POSITIVE_INFINITY;
            PruningControlSuper pruning3 = cfp.setupPruning(searchSpace, pruningInterval3, false, false);
            pruning3.prune();//pass in DEE options, and run the specified types of DEE 

            SelfConsistentMeanField_Parallel scmf2 = new SelfConsistentMeanField_Parallel(mrf);
            scmf2.run();
            BigDecimal ZLB2 = scmf2.calcPartitionFunction();
            double logZLB_2 = ef.log10(ZLB2);
            System.out.println("Lower bound on log partition function (SCMF Parallel) = " + logZLB_2);

            double numConfs = 1;
            for (int pos = 0; pos < searchSpace.emat.numPos(); pos++) {
                numConfs = numConfs * searchSpace.pruneMat.unprunedRCsAtPos(pos).size();
            }
            //BigDecimal Zpart = calcRigidPartFunction(searchSpace);
            //BigDecimal logZpart = ef.log(Zpart);
            MapPerturbation mapPert = new MapPerturbation(searchSpace);
            double logZUB = (0.4342944819) * mapPert.calcUBLogZ(100);
            System.out.println("Upper bound on log partition function (MAP-Pert) = " + logZUB);
            ArrayList<Integer> toMerge = mapPert.getPairWithMaxMutualInfo(true);
            searchSpace.mergePositionRigid(toMerge);
            searchSpace.competitorPruneMat = null;

            try (PrintStream out = new PrintStream(new FileOutputStream("results.txt", true))) {
                out.println("GMEC: " + Double.toString(gmecScore));
                out.println("Lower Bound SCMF: " + Double.toString(logZLB));
                out.println("Lower Bound SCMF_Parallel: " + Double.toString(logZLB_2));
                out.println("Upper Bound: " + Double.toString(logZUB));
            } catch (Exception e) {
            }
        }
        /*
         try( PrintStream out = new PrintStream(new FileOutputStream("results.txt", false)) ) {
         out.println("GMEC: "+Double.toString(gmecE));
         out.println("Lower Bound SCMF: "+Double.toString(logZLB));
         out.println("Lower Bound SCMF_Parallel: "+Double.toString(logZLB_2));
         out.println("Upper Bound: "+Double.toString(logZUB));
         }
         catch(Exception e){}
         */

        /*
         ArrayList<Integer> posToMerge = new ArrayList<>();
         posToMerge.add(3);
         posToMerge.add(5);
         searchSpace.mergePositionContinuous(posToMerge);
         if (searchSpace.competitorPruneMat == null) {
         System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
         PruningControlSuper compPruning = cfp.setupPruning(searchSpace, 0, false, false);
         compPruning.setOnlyGoldstein(true);
         compPruning.prune();
         searchSpace.competitorPruneMat = searchSpace.pruneMat;
         searchSpace.pruneMat = null;
         System.out.println("COMPETITOR PRUNING DONE");
         }
         double pruningInterval2 = 30.0;
         PruningControlSuper pruning2 = cfp.setupPruning(searchSpace, pruningInterval2, false, false);
         pruning2.prune();//pass in DEE options, and run the specified types of DEE 
         double gmecE2 = -(0.4342944819) * calcGMEC(searchSpace) / this.constRT;
         */
    }

    //getGMEC from lower bounds
    private BigDecimal calcRigidPartFunction(SearchProblemSuper aSearchSpace) {
        boolean needToRepeat;
        int[] GMECConf = null;
        double bestESoFar = Double.POSITIVE_INFINITY;

        SearchProblemSuper searchSpace = aSearchSpace;
        BigDecimal partFunction = new BigDecimal(0.0);

        int iter = 0;
        ConfSearch search = new ConfTreeSuper(searchSpace);

        do {
            needToRepeat = true;
            int[] conf = search.nextConf();
            double E = searchSpace.lowerBound(conf);
            partFunction = partFunction.add(ef.exp(-(E) / constRT));
            iter++;
            if (iter > 0) {
                needToRepeat = false;
            }
        } while (needToRepeat);
        return partFunction;
    }

    //getGMEC from lower bounds
    private double calcGMEC(SearchProblemSuper aSearchSpace) {
        SearchProblemSuper searchSpace = aSearchSpace;
        ConfSearch search = new ConfTreeSuper(searchSpace);
        int[] conf = search.nextConf();
        double E = searchSpace.lowerBound(conf);
        return E;
    }

}
