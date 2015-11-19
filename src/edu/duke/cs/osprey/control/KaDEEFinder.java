/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.util.List;

import edu.duke.cs.osprey.astar.ConfTreeSuper;
import edu.duke.cs.osprey.astar.comets.*;
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
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.StringParsing;
import java.awt.Point;
import java.io.PrintStream;
import java.io.FileOutputStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.stream.Collectors;
import java.util.List;
import java.util.ArrayList;
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
    //0: Bound 
    //1: Mutable UnBound
    //2: Mutable Bound
    SearchProblemSuper[] searchSpaces;

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
        this.searchSpaces = cfp.getSearchProblemSupers();

        if (true) {
            COMETSTreeSuper tree = setupCometsTree();
            int[] seq1 = tree.nextConf();
            //while (true){
            //    int[] seq2 = tree.nextConf();
           // }
        } else {

            SearchProblemSuper searchSpace = searchSpaces[0];
            searchSpace.loadEnergyMatrix();
            //Doing competitor pruning now
            //will limit us to a smaller, but effective, set of competitors in all future DEE
            List<Integer> mut2StatePosNumBound = new ArrayList<>();
            List<Integer> mut2StatePosNumUnBound = new ArrayList<>();
            int numMutRes = 0;

            ConfSpaceSuper cSpaceMonomerMut = searchSpaces[1].confSpaceSuper;
            ArrayList<ArrayList<String>> allowedAAComplex = cfp.getAllowedAAs();
            for (int i = 0; i < cSpaceMonomerMut.posFlexSuper.size(); i++) {
                if (allowedAAComplex.get(i).size() > 0) {
                    mut2StatePosNumBound.add(i);
                    mut2StatePosNumUnBound.add(i);
                    numMutRes++;
                }
            }
            List<List<Integer>> mutable2StatePosNums = new ArrayList<>();
            mutable2StatePosNums.add(mut2StatePosNumBound);
            mutable2StatePosNums.add(mut2StatePosNumUnBound);

            LME objFcn = new LME(new double[]{1, -1}, 0.0, 2);

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
            double pruningInterval = 0.0;
            PruningControlSuper pruning = cfp.setupPruning(searchSpace, pruningInterval, false, false);
            pruning.prune();//pass in DEE options, and run the specified types of DEE 
//        BigInteger confSpace = new BigInteger("1");
//        for (int pos = 0; pos < searchSpace.emat.numPos(); pos++){
//            confSpace = confSpace.multiply(new BigInteger(((Integer) searchSpace.emat.oneBody.get(pos).size()).toString()));
//        }
            //System.out.println(confSpace);

            double gmecE = calcGMEC(searchSpace);
            System.out.println("GMEC Energy = " + gmecE);
            //Calculate GMEC (0.43429 = log_10(e) puts gmecE in log_10 scale 
            double gmecScore = -(0.4342944819) * gmecE / this.constRT;
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
            List<Integer> toMerge = mapPert.getPairWithMaxMutualInfo(true);
            logZLB = (0.4342944819) * mapPert.calcLBLogZ(100);
            System.out.println("Lower bound on log partition function (MAP-Pert) = " + logZLB);

            //searchSpace.mergePositionRigid(toMerge);
            //searchSpace.competitorPruneMat = null;
            //searchSpace.pruneMat = null;
            try (PrintStream out = new PrintStream(new FileOutputStream("results.txt", true))) {
                out.println("GMEC: " + Double.toString(gmecScore));
                out.println("Lower Bound SCMF: " + Double.toString(logZLB));
                out.println("Lower Bound SCMF_Parallel: " + Double.toString(logZLB_2));
                out.println("Upper Bound: " + Double.toString(logZUB));
            } catch (Exception e) {
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
             List<Integer> posToMerge = new ArrayList<>();
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

    //Loads energy matrices and prune for COMETS
    private void loadEMatandPruneComets(double pruningInterval) {
        cfp.params.setValue("TYPEDEP", "TRUE");
        cfp.params.setValue("BOUNDSTHRESH", "100000000000000");
        for (int state = 0; state < searchSpaces.length; state++) {
            SearchProblemSuper searchProblem = this.searchSpaces[state];

            System.out.println("Precomputing Energy Matrix for " + searchProblem.name + " state");
            searchProblem.loadEnergyMatrix();

            System.out.println("Initializing Pruning for " + searchProblem.name + " state");
            initializePruning(searchProblem);
            PruningControlSuper pruning = cfp.setupPruning(searchProblem, pruningInterval, useEPIC, useTupExp);
            pruning.prune();
        }
    }

    //Given three search problems (Bound, UnBound Prot, Unbound Lig) this function
    //sets up the COMETS tree.
    //The nonmutable unbound state is added and used just as a constant to the objective function
    private COMETSTreeSuper setupCometsTree() {

        //For each state, for each position, this contains a list of allowed 
        //amino acids at that position
        List<List<List<String>>> allowedAAsPerState = new ArrayList<>();
        for (int state = 0; state < searchSpaces.length; state++) {
            allowedAAsPerState.add(getAllowedAA(state));
        }

        //Load the energy matrices and do pruning
        double pruningInterval = 0.0;
        loadEMatandPruneComets(pruningInterval);
        /*
         for (SearchProblemSuper searchProblem : searchSpaces){
         searchProblem.emat.setConstTerm(0.0);
         }
         */
        //For each state this arraylist gives the mutable pos nums of that state
        List<List<Integer>> mutable2StatePosNum = handleMutable2StatePosNums(allowedAAsPerState);

        //determine which states are mutable and which are non-mutable
        boolean[] stateIsMutable = new boolean[this.searchSpaces.length];
        int numMutableState = 0;
        int numNonMutableState = 0;
        for (int state = 0; state < searchSpaces.length; state++) {
            stateIsMutable[state] = !(mutable2StatePosNum.get(state).isEmpty()); //If not empty, it is mutable
            if (stateIsMutable[state]) {
                numMutableState++;
            } else {
                numNonMutableState++;
            }
        }

        //Get Mutable States
        SearchProblemSuper[] mutableStates = new SearchProblemSuper[numMutableState];
        //For each MUTABLE state, this arraylist gives the mutable pos nums of that state
        //This is what will go into the COMETS tree
        List<List<Integer>> mutableState2StatePosNumList = new ArrayList<>();
        List<List<List<String>>> mutableStateAllowedAAs = new ArrayList<>();
        int mutableStateIndex = 0;

        SearchProblemSuper nonMutableState = searchSpaces[1];
        double unboundLigandGMECEnergy = 0.0;
        for (int state = 0; state < searchSpaces.length; state++) {
            if (stateIsMutable[state]) {
                mutableStates[mutableStateIndex] = searchSpaces[state];
                mutableState2StatePosNumList.add(mutable2StatePosNum.get(state));
                mutableStateAllowedAAs.add(allowedAAsPerState.get(state));
                mutableStateIndex++;
            } else {
                nonMutableState = searchSpaces[state];
                //For the const term of the LME objective function
                unboundLigandGMECEnergy = getGMECEnergyLigand(nonMutableState);
            }
        }
        //For const term of LME objective function
        int numStatesForCOMETS = mutableStates.length;
        int numTreeLevels = getNumMutablePos(mutableState2StatePosNumList);
        ArrayList<ArrayList<String>> AATypeOptions = handleAATypeOptions(mutableStateAllowedAAs);
        LME objFcn = new LME(new double[]{1, -1}, -unboundLigandGMECEnergy, 2);
        LME[] constraints = new LME[0];
        int numMaxMut = -1;
        String[] wtSeq = null;

        //Convert to ArrayList<>
        ArrayList<ArrayList<Integer>> mutableState2StatePosNum = new ArrayList<>();
        for (List<Integer> mutable2PosNum : mutableState2StatePosNumList) {
            ArrayList<Integer> converted = new ArrayList(mutable2PosNum);
            mutableState2StatePosNum.add(converted);
        }
        int[] numMutPerStrand = cfp.getNumMutPerStrand();

        COMETSTreeSuper tree = new COMETSTreeSuper(numTreeLevels, objFcn, constraints, AATypeOptions, numMaxMut, wtSeq, mutableStateIndex, mutableStates, nonMutableState, mutableState2StatePosNum, numMutPerStrand);
        return tree;
    }

    private List<List<String>> getAllowedAA(int state) {
        ArrayList<ArrayList<String>> complexAllowedAAs = cfp.getAllowedAAs();
        ArrayList<String> complexFlexRes = cfp.getFlexRes();
        int numFlexRes = complexFlexRes.size();
        int beginPos = -1;
        int endPos = -1;
        if (state == 0) {
            beginPos = 0;
            endPos = numFlexRes;
        } else if (state == 1) {
            beginPos = 0;
            endPos = this.searchSpaces[1].confSpaceSuper.posFlexSuper.size();
        } else {
            beginPos = this.searchSpaces[1].confSpaceSuper.posFlexSuper.size();
            endPos = beginPos + this.searchSpaces[2].confSpaceSuper.posFlexSuper.size();
        }
        List<List<String>> allowedAAs = new ArrayList<>();
        Molecule wtMolec = PDBFileReader.readPDBFile(cfp.params.getValue("PDBNAME"));

        for (int posNum = beginPos; posNum < endPos; posNum++) {
            List<String> currentAAOptions = complexAllowedAAs.get(posNum);
            if (cfp.params.getBool("AddWT", true)) {
                Residue res = wtMolec.getResByPDBResNumber(complexFlexRes.get(posNum));
                if (!StringParsing.containsIgnoreCase(complexFlexRes, res.template.name)) {
                    currentAAOptions.add(res.template.name);
                }
            }
            if (currentAAOptions.isEmpty()) {
                throw new RuntimeException("ERROR: No AAtype for Position: " + posNum);
            }
            allowedAAs.add(currentAAOptions);
        }
        return allowedAAs;
    }

    private String[] getWTSequence(List<Integer> complexMut2StatePosNum) {
        String[] wtSequence = new String[complexMut2StatePosNum.size()];
        ArrayList<String> complexFlexRes = cfp.getFlexRes();
        Molecule wtMolec = PDBFileReader.readPDBFile(cfp.params.getValue("PDBNAME"));

        for (int posIndex = 0; posIndex < complexMut2StatePosNum.size(); posIndex++) {
            Integer mutPosNum = complexMut2StatePosNum.get(posIndex);
            String mutPos = complexFlexRes.get(mutPosNum);
            wtSequence[posIndex] = wtMolec.getResByPDBResNumber(mutPos).template.name;
        }
        return wtSequence;
    }

    //for each state, returns a list of mutable pos numbers. 
    //If the state has no mutable positions an empty list is returned
    //Thus, we can use this to decide which states go in the comets search
    private List<List<Integer>> handleMutable2StatePosNums(List<List<List<String>>> allowedAAPerState) {
        List<List<Integer>> mutable2StatePosNum = new ArrayList<>();

        for (int state = 0; state < allowedAAPerState.size(); state++) {
            List<List<String>> allowedAAForState = allowedAAPerState.get(state);
            List<Integer> mutablePositionsForState = getMutablePosNums(allowedAAForState);
            mutable2StatePosNum.add(mutablePositionsForState);
        }

        return mutable2StatePosNum;
    }

    //Given the allowed AAs at a state, get the mutable position numbers
    private List<Integer> getMutablePosNums(List<List<String>> allowedAAForState) {
        List<Integer> mutablePosNum = new ArrayList<>();
        for (int posNum = 0; posNum < allowedAAForState.size(); posNum++) {
            if (allowedAAForState.get(posNum).size() > 1) {
                mutablePosNum.add(posNum);
            }
        }
        return mutablePosNum;
    }

    private int getNumMutablePos(List<List<Integer>> mutable2StatePosNum) {
        int numMutable = -1;
        for (int state = 0; state < mutable2StatePosNum.size(); state++) {
            int numMutableAtState = mutable2StatePosNum.get(state).size();
            numMutable = Math.max(numMutable, numMutableAtState);
        }
        return numMutable;
    }

    //Make sure each mutable state has the same Allowed AA
    //Return Allowed AA for all mutable positions
    private ArrayList<ArrayList<String>> handleAATypeOptions(List<List<List<String>>> mutableStateAllowedAAs) {
        List<List<String>> AATypeOptionsList = mutableStateAllowedAAs.get(0).stream().filter(aaTypes -> aaTypes.size() > 1).collect(Collectors.toList());
        //Convert to ArrayList<>
        ArrayList<ArrayList<String>> AATypeOptions = new ArrayList<>();
        for (List<String> AAOption : AATypeOptionsList) {
            ArrayList<String> converted = new ArrayList<String>(AAOption);
            AATypeOptions.add(converted);
        }

        for (int state = 1; state < mutableStateAllowedAAs.size(); state++) {
            List<List<String>> AATypesForState = mutableStateAllowedAAs.get(state);
            if (AATypesForState.size() != AATypeOptions.size()) {
                throw new RuntimeException("ERROR: Different Number of Mutable Positions between Bound and Unbound");
            }
            for (int posNum = 0; posNum < AATypesForState.size(); posNum++) {
                List<String> AATypesForPos = AATypesForState.get(posNum);
                for (int aaIndex = 0; aaIndex < AATypesForPos.size(); aaIndex++) {
                    if (!(AATypeOptions.get(posNum).contains(AATypesForPos.get(aaIndex)))) {
                        throw new RuntimeException("ERROR: AAType Different for Bound and Unbound Mutable Residues");
                    }
                }
            }
        }
        return AATypeOptions;
    }

    private double getGMECEnergyLigand(SearchProblemSuper searchProblem) {
        ConfTreeSuper confTree = new ConfTreeSuper(searchProblem);
        int[] conf = confTree.nextConf();
        SuperRCTuple rcTup = new SuperRCTuple(conf);
        double unboundLigandE = searchProblem.emat.getInternalEnergy(rcTup) + searchProblem.emat.getConstTerm();
        return unboundLigandE;
    }

    private void initializePruning(SearchProblemSuper searchProblem) {
        //Creates an efficient competitor pruning matrix
        searchProblem.competitorPruneMat = null;
        System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
        //prune with 0 interval, anything that survives will be added as a competitor
        PruningControlSuper compPruning = cfp.setupPruning(searchProblem, 0, false, false);
        compPruning.setOnlyGoldstein(true);
        compPruning.prune();
        searchProblem.competitorPruneMat = searchProblem.pruneMat;
        searchProblem.pruneMat = null;
        System.out.println("COMPETITOR PRUNING DONE");
    }
}
