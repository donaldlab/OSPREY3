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
import edu.duke.cs.osprey.markovrandomfield.MarkovRandomField;
import edu.duke.cs.osprey.markovrandomfield.SelfConsistentMeanField;
import edu.duke.cs.osprey.markovrandomfield.GumbelDistribution;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import java.awt.Point;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.stream.Collectors;
import java.util.List;
import javafx.scene.input.KeyCode;

/**
 *
 * @author hmn5
 */
public class KaDEEFinder {

    ConfigFileParser cfp;

    SearchProblemSuper searchSpace;

    double Ew;//energy window for enumerating conformations: 0 for just GMEC
    double I0 = 0;//initial value of iMinDEE pruning interval
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

    void doKaDEE() {
        //Calculate the GMEC
        ArrayList<Integer> testing = new ArrayList<>();
        testing.add(1);
        testing.add(2);
        testing.add(3);
        List<Integer> h = testing.stream().map(x -> x + 1).collect(Collectors.toList());

        double curInterval = I0;//For iMinDEE.  curInterval will need to be an upper bound

        searchSpace = cfp.getSearchProblemSuper();
        ConfSpaceSuper confSpaceSuper = searchSpace.confSpaceSuper;
        searchSpace.loadEnergyMatrix();

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

        MarkovRandomField mrf = new MarkovRandomField(searchSpace, 0.0);
        SelfConsistentMeanField scmf = new SelfConsistentMeanField(mrf);
        scmf.run();
        BigDecimal Z = scmf.calcPartitionFunction();
        BigDecimal logZ = ef.log(Z);

        double numConfs = 1;
        for (int pos = 0; pos < searchSpace.emat.numPos(); pos++) {
            numConfs = numConfs * searchSpace.pruneMat.unprunedRCsAtPos(pos).size();
        }
        //BigDecimal Zpart = calcRigidPartFunction(searchSpace);
        //BigDecimal logZpart = ef.log(Zpart);
        double logZUB = calcLogPartMapPert(searchSpace);
        searchSpace.loadEnergyMatrix();
        BigDecimal Zpart = calcRigidPartFunction(searchSpace);
        searchSpace.loadEnergyMatrix();
        
        
        
        ArrayList<ArrayList<Integer>> posToMerge = new ArrayList<>();
        for (int i = 0; i < confSpaceSuper.posFlex.size(); i++) {
            ArrayList<Integer> newPos = new ArrayList<>();
            if (i == 0) {
                newPos.add(i);
                i++;
                newPos.add(i);
            } else {
                newPos.add(i);
            }
            posToMerge.add(newPos);
        }
        searchSpace.mergePositionRigid(posToMerge.get(0));

        double pruningInterval2 = Double.POSITIVE_INFINITY;
        //Doing competitor pruning now
        //will limit us to a smaller, but effective, set of competitors in all future DEE
        searchSpace.competitorPruneMat = null;
        searchSpace.pruneMat = null;
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
        PruningControlSuper pruning2 = cfp.setupPruning(searchSpace, pruningInterval2, false, false);
        pruning2.prune();//pass in DEE options, and run the specified types of DEE 

        BigDecimal ZpartMerges = calcRigidPartFunction(searchSpace);
        double logZUBMerged = calcLogPartMapPert(searchSpace);

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
    private double calcLogPartMapPert(SearchProblemSuper aSearchSpace) {
        boolean needToRepeat;
        int[] GMECConf = null;
        double bestESoFar = Double.POSITIVE_INFINITY;
        int numSamples = 100;
        SearchProblemSuper searchSpace = aSearchSpace;
        BigDecimal averageGMECs = new BigDecimal(0.0);
        int iter = 0;
        for (int i = 0; i < numSamples; i++) {
            //searchSpace.loadMergedEnergyMatrix();
            ArrayList<ArrayList<Double>> originalOneBodyEmat = (ArrayList<ArrayList<Double>>) ObjectIO.deepCopy(this.searchSpace.emat.oneBody);
            addGumbelNoiseOneBody(searchSpace);
            ConfSearch search = new ConfTreeSuper(searchSpace);
            int[] conf = search.nextConf();
            double E = -1.0 * searchSpace.lowerBound(conf);
            averageGMECs = averageGMECs.add(new BigDecimal(E));
            //replace oneBody with original to remove the noise added
            this.searchSpace.emat.oneBody = originalOneBodyEmat;
            iter++;
        }
        return averageGMECs.divide(new BigDecimal(numSamples * this.constRT),ef.mc).doubleValue();
    }

    //add Gumbel noise to one-body terms
    private void addGumbelNoiseOneBody(SearchProblemSuper aSearchSpace) {
        SearchProblemSuper searchSpace = aSearchSpace;
        EnergyMatrix emat = searchSpace.emat;
        for (int pos = 0; pos < emat.oneBody.size(); pos++) {
            for (int superRC : searchSpace.pruneMat.unprunedRCsAtPos(pos)) {
                double currentE = emat.getOneBody(pos, superRC);
                double noise = GumbelDistribution.sample(-1.0*GumbelDistribution.gamma,1.0) * this.constRT;
                emat.setOneBody(pos, superRC, currentE - noise);
            }
        }
    }
    
}
