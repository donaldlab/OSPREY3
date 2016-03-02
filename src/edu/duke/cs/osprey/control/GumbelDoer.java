/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;


import edu.duke.cs.osprey.astar.kadee.GumbelMapTree;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.partitionfunctionbounds.MapPerurbation;
import edu.duke.cs.osprey.partitionfunctionbounds.MarkovRandomField;
import edu.duke.cs.osprey.partitionfunctionbounds.SelfConsistentMeanField;
import edu.duke.cs.osprey.pruning.PruningControl;

/**
 *
 * @author hmn5
 */
public class GumbelDoer {
    ConfigFileParser cfp;
    SearchProblem sp;
    final double constRT = PoissonBoltzmannEnergy.constRT;

    public GumbelDoer(ConfigFileParser aCFP) {
        this.cfp = aCFP;
        sp = cfp.getSearchProblem();
        loadEMatandPrune(Double.POSITIVE_INFINITY);
        
        double average = 0.0;
        int averageNodesExpanded = 0;
        int numSamples = 50;
        for (int i = 0; i < numSamples; i++) {
            GumbelMapTree tree = new GumbelMapTree(sp);
            tree.nextConf();
            double score = tree.currentBestFeasibleScore;
            average += score;
            double logZestimate = -average/(this.constRT*(i+1));
            System.out.println("Number of Nodes Expanded: "+tree.numExpanded);
            averageNodesExpanded += tree.numExpanded;
            System.out.println("Current Sample: "+-score/this.constRT);
            System.out.println("Current Average: "+logZestimate);
            System.out.println("Current Average Nodes Exp: "+averageNodesExpanded/(i+1));
        }
        System.out.println("Average Nodes Expanded: "+averageNodesExpanded/numSamples);
         double logZ = -average / (this.constRT * numSamples);
         System.out.println("Gumbel logZ: "+logZ);
         
        MarkovRandomField mrf = new MarkovRandomField(sp, 0.0);
        SelfConsistentMeanField scmf = new SelfConsistentMeanField(mrf);
        scmf.run();
        double lowerZ = scmf.calcLBLogZ();
        System.out.println("SCMF LogZ: "+lowerZ);
        MapPerurbation mpert = new MapPerurbation(sp);
        double upperZ = mpert.calcUBLogZ(500);
        System.out.println("MapPert LogZ:" + upperZ);
    }

    //Loads energy matrices and prune 
    private void loadEMatandPrune(double pruningInterval) {
        System.out.println("Precomputing Energy Matrix for " + sp.name + " state");
        sp.loadEnergyMatrix();

        System.out.println("Initializing Pruning for " + sp.name + " state");
        initializePruning(sp);
        PruningControl pruning = cfp.setupPruning(sp, pruningInterval, false, false);
        pruning.prune();
    }

    private void initializePruning(SearchProblem searchProblem) {
        //Creates an efficient competitor pruning matrix
        searchProblem.competitorPruneMat = null;
        System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
        //prune with 0 interval, anything that survives will be added as a competitor
        PruningControl compPruning = cfp.setupPruning(searchProblem, Double.POSITIVE_INFINITY, false, false);
        compPruning.setOnlyGoldstein(true);
        compPruning.prune();
        searchProblem.competitorPruneMat = searchProblem.pruneMat;
        searchProblem.pruneMat = null;
        System.out.println("COMPETITOR PRUNING DONE");
    }
}

