/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.TermECalculatorSuper;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.ConfETupleExpander;
import edu.duke.cs.osprey.tupexp.TupExpChooser;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.stream.Collectors;

/**
 *
 * @author hmn5
 */
public class SearchProblemSuper implements Serializable{
    //This is the super-RC analog of SearchProblem

    public ConfSpaceSuper confSpaceSuper;
    public ConfSpace confSpace;

    public EnergyMatrix emat;

    public EPICMatrix epicMat = null;
    public EPICSettings epicSettings = null;

    public EnergyMatrix tupExpEMat;

    public EnergyFunction fullConfE;//full energy for any conformation
    public ArrayList<Residue> shellResidues;

    public String name;

    public PruningMatrix pruneMat;

    boolean contSCFlex;

    public PruningMatrix competitorPruneMat;

    public boolean useEPIC = false;
    public boolean useTupExpForSearch = false;

    public SearchProblemSuper(SearchProblemSuper sp1) {
        //Shallow copy

        confSpaceSuper = sp1.confSpaceSuper;
        emat = sp1.emat;
        epicMat = sp1.epicMat;
        epicSettings = sp1.epicSettings;
        tupExpEMat = sp1.tupExpEMat;

        fullConfE = sp1.fullConfE;
        shellResidues = sp1.shellResidues;
        name = sp1.name + System.currentTimeMillis();//probably will want to change this to something more meaningful

        pruneMat = sp1.pruneMat;
        competitorPruneMat = sp1.competitorPruneMat;

        contSCFlex = sp1.contSCFlex;
        useEPIC = sp1.useEPIC;
        useTupExpForSearch = sp1.useTupExpForSearch;
    }

    public SearchProblemSuper(String name, String PDBFile, ArrayList<String> flexibleRes, ArrayList<ArrayList<String>> allowedAAs, boolean addWT,
            boolean addWTRots, boolean contSCFlex, boolean useEPIC, EPICSettings epicSettings, boolean useTupExp, DEEPerSettings dset,
            ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, boolean useEllipses) {

        this.confSpaceSuper = new ConfSpaceSuper(PDBFile, flexibleRes, allowedAAs, addWT, addWTRots, contSCFlex, dset, moveableStrands, freeBBZones, useEllipses);
        this.confSpace = confSpaceSuper;
        this.name = name;

        this.contSCFlex = contSCFlex;
        this.useTupExpForSearch = useTupExp;
        this.useEPIC = useEPIC;
        this.epicSettings = epicSettings;

        //energy function setup
        EnergyFunctionGenerator eGen = EnvironmentVars.curEFcnGenerator;
        decideShellResidues(eGen.distCutoff);
        //DEBUG
        //addAllResToShell();
        //fullConfE = eGen.fullConfEnergy(confSpace,shellResidues);
        //HMN: use confSpaceSuper as input instead of confSpace
        fullConfE = eGen.fullConfEnergy(confSpaceSuper, shellResidues);
    }

    public SearchProblemSuper getSubsetSearchProblem(ArrayList<Integer> subsetOfPositions) {
        Collections.sort(subsetOfPositions);

        SearchProblemSuper subsetSP = new SearchProblemSuper(this);
        if (useEPIC || useTupExpForSearch) {
            throw new RuntimeException("ERROR: Currently getSubsetSearchProblem does not support EPIC or LUTE");
        }
        if (contSCFlex) {
            //we will want to create a new energy function generator that properly
            //handles this partial space 
        }

        
        subsetSP.confSpaceSuper = this.confSpaceSuper.getSubsetConfSpace(subsetOfPositions);
        
        // Begin: Create new Energy and Pruning Matrix (currently only one and two-body terms)
        EnergyMatrix subsetEmat = new EnergyMatrix(subsetSP.confSpaceSuper, this.emat.getPruningInterval());
        ArrayList<ArrayList<Double>> newOneBodyE = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> newPairwiseE = new ArrayList<>();

        PruningMatrix subsetPruneMat = new PruningMatrix(confSpaceSuper, this.pruneMat.getPruningInterval());
        ArrayList<ArrayList<Boolean>> newOneBodyPrune = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<ArrayList<Boolean>>>> newPairwisePrune = new ArrayList<>();

        for (int i = 0; i < subsetSP.confSpaceSuper.numPos; i++) {
            int posNumI = subsetOfPositions.get(i);
            ArrayList<Double> oneBodyEAtPosI = (ArrayList<Double>) ObjectIO.deepCopy(this.emat.oneBody.get(posNumI));
            ArrayList<Boolean> oneBodyPruneAtPosI = (ArrayList<Boolean>) ObjectIO.deepCopy(this.pruneMat.oneBody.get(posNumI));
            ArrayList<ArrayList<ArrayList<Double>>> pairwiseEAtPosI = new ArrayList<>();
            ArrayList<ArrayList<ArrayList<Boolean>>> pairwisePruneAtPosI = new ArrayList<>();
            for (int j = 0; j < i; j++) {
                int posNumJ = subsetOfPositions.get(j);
                ArrayList<ArrayList<Double>> pairwiseEAtPosIPosJ = (ArrayList<ArrayList<Double>>) ObjectIO.deepCopy(this.emat.pairwise.get(posNumI).get(posNumJ));
                pairwiseEAtPosI.add(pairwiseEAtPosIPosJ);
                ArrayList<ArrayList<Boolean>> pairwisePruneAtPosIPosJ = (ArrayList<ArrayList<Boolean>>) ObjectIO.deepCopy(this.pruneMat.pairwise.get(posNumI).get(posNumJ));
                pairwisePruneAtPosI.add(pairwisePruneAtPosIPosJ);
            }
            newOneBodyE.add(oneBodyEAtPosI);
            newPairwiseE.add(pairwiseEAtPosI);
            newOneBodyPrune.add(oneBodyPruneAtPosI);
            newPairwisePrune.add(pairwisePruneAtPosI);
        }
        subsetEmat.oneBody = newOneBodyE;
        subsetEmat.pairwise = newPairwiseE;
        //For now lets not handle higher-order terms
        
        subsetSP.emat = subsetEmat;
        subsetPruneMat.oneBody = newOneBodyPrune;
        subsetPruneMat.pairwise = newPairwisePrune;
        subsetSP.pruneMat = subsetPruneMat;
        // End: Create new Energy and Pairwsie Matrix
        
        return subsetSP;
    }

    private void addAllResToShell() {
        shellResidues = new ArrayList<>();
        for (Residue nonFlexRes : confSpaceSuper.m.residues) {
            shellResidues.add(nonFlexRes);
        }
    }

    private void decideShellResidues(double distCutoff) {
        //Decide what non-flexible residues need to be accounted for in energy calculations

        //use array list for these so they stay in order
        ArrayList<Residue> flexibleResidues = new ArrayList<>();
        for (PositionConfSpaceSuper pcs : confSpaceSuper.posFlexSuper) {
            for (Residue res : pcs.resList) {
                flexibleResidues.add(res);
            }
        }

        //we'll decide what shell residues to include by using a simple distance cutoff with
        //the current conformation,
        //rather than doing a conformational search for the minimal distance (with respect to conformations)
        //of a shell residue to any flexible residues
        //the distance cutoff can be increased to accommodate this if desired.
        shellResidues = new ArrayList<>();

        for (Residue nonFlexRes : confSpaceSuper.m.residues) {
            //if residue is not flexible
            if (!flexibleResidues.contains(nonFlexRes)) {
                for (Residue flexRes : flexibleResidues) {
                    double dist = flexRes.distanceTo(nonFlexRes);
                    if (dist <= distCutoff) {
                        shellResidues.add(nonFlexRes);//close enough to flexRes that we should add it
                        break;
                    }
                }
            }
        }
    }

    //We want to make sure the monomer has the same shell residues as the complex
    //that are also in the strand
    //Given: Shell Residues from Complex
    //Match Residues in this molecule with those in the complex molecule 
    /////(note the residues are different objects so we compare using fullname)
    //Output: Shell residues from this molecule
    public void getShellResiduesFromComplex(ArrayList<Residue> shellResFromComplex) {
        ArrayList<Residue> newShellResidues = new ArrayList<>();

        //Get flexible residue
        ArrayList<Residue> flexibleResidues = new ArrayList<>();
        for (PositionConfSpaceSuper pcs : confSpaceSuper.posFlexSuper) {
            for (Residue res : pcs.resList) {
                flexibleResidues.add(res);
            }
        }

        for (Residue res : this.confSpaceSuper.m.residues) {
            if (!flexibleResidues.contains(res)) {
                for (Residue complexShellRes : shellResFromComplex) {
                    if (complexShellRes.fullName.equals(res.fullName)) {
                        newShellResidues.add(res);
                        break;
                    }
                }
            }
        }
        //Are the new shell residues different from what we originall had
        boolean newShellIsDifferent = false;
        for (Residue shellRes : newShellResidues) {
            if (!this.shellResidues.contains(shellRes)) {
                newShellIsDifferent = true;
            }
        }
        for (Residue oldShellRes : this.shellResidues) {
            if (!newShellResidues.contains(oldShellRes)) {
                newShellIsDifferent = true;
            }
        }
        //if the shell residues are different calculate new fullConfE
        if (newShellIsDifferent) {
            this.shellResidues = newShellResidues;
            EnergyFunctionGenerator eGen = EnvironmentVars.curEFcnGenerator;
            this.fullConfE = eGen.fullConfEnergy(this.confSpaceSuper, shellResidues);
        }
    }

    public double minimizedEnergy(int[] conf) {
        //Minimized eneryg of the conformatoin
        //whose super-RCs are listed for all flexible positions in conf
        double E = confSpaceSuper.minimizeEnergy(conf, fullConfE, null);
        //HMN: Add constant term
        E += emat.getConstTerm();
        return E;
    }

    public void outputMinimizedStruct(int[] conf, String PDBFileName) {
        //Output minimized conformation to specified file
        //super-RCs are listed for all flexible positions in conf
        confSpaceSuper.minimizeEnergy(conf, fullConfE, PDBFileName);
    }

    //HMN: Not currently supported
    public double approxMinimizedEnergy(int[] conf) {
        //EPIC or other approximation for the minimized energy of the conformation
        //whose RCs are listed for all flexible positions in conf
        throw new RuntimeException("TupExpansion and EPIC not currently supported with SearchProblemSuper and super-RCs");
        /*
         if (useTupExpForSearch) {//use tup-exp E-matrix direectly
         return tupExpEMat.confE(conf);
         } else if (useEPIC) {//EPIC w/o tup-exp
         return EPICMinimizedEnergy(conf);
         } else {
         throw new RuntimeException("ERROR: Asking searchSpace to approximate minimized energy but using a non-approximating method");
         }
         */
    }

    //HMN: Not currently supported
    public double EPICMinimizedEnergy(int[] conf) {
        throw new RuntimeException("TupExpansion and EPIC not currently supported with SearchProblemSuper and super-RCs");
        /*
         //approximate energy using EPIC
         double bound = emat.confE(conf);//emat contains the pairwise lower bounds
         double contPart = epicMat.minContE(conf);
         //EPIC handles the continuous part (energy - pairwise lower bounds)

         return bound + contPart;
         */
    }

    public double lowerBound(int[] conf) {
        //Argument: RC assignments for all the flexible residues (RCs defined in resFlex)
        //return lower bound on energy for the conformational space defined by these assignments
        //based on precomputed energy matrix (including EPIC if indicated)

        double bound = emat.confE(conf);//the energy recorded by the matrix is 
        //the pairwise lower bounds

        return bound;
    }

    //LOADING AND PRECOMPUTATION OF ENERGY MATRIX-TYPE OBJECTS (regular energy matrix, tup-exp and EPIC matrices)
    public void loadEnergyMatrix() {
        loadMatrix(MatrixType.EMAT);

    }

    public void loadTupExpEMatrix() {
        loadMatrix(MatrixType.TUPEXPEMAT);
    }

    public void loadEPICMatrix() {
        loadMatrix(MatrixType.EPICMAT);
    }

    enum MatrixType {

        EMAT, TUPEXPEMAT, EPICMAT;
    }

    //This function will merge the positions and then merge the energy matrix,
    //recomputing only those terms that need to be recomputed
    //posToMerge contains all the positions to merges
    public void mergePositionContinuous(ArrayList<Integer> posToMerge) {
        //newPosList contains all the info about our new merged position
        //In this example newPosList will be [[0],[1,3],[2],[4]]
        ArrayList<ArrayList<Integer>> newPosList = getNewPosList(posToMerge);
        //get the index in the list that has multiple positions
        int newPosIndex = getNewPosIndex(newPosList);

        //If the number of positions to merge is 2, we can merge
        if (posToMerge.size() == 2) {
            //Create temp prune Mat to decide which RCs to merge
            PruningMatrix newPruneMat = mergePruneMat(newPosList);
            confSpaceSuper.mergePosition(newPosList, newPruneMat);
            this.mergeContinuousEnergyMatrix(newPosList, newPruneMat);
        } //Otherwise, we recurse by merging the first two positoin
        else {
            ArrayList<Integer> toMergeNow = new ArrayList<>();
            toMergeNow.add(posToMerge.get(0));
            toMergeNow.add(posToMerge.get(1));
            ArrayList<Integer> remainingPos = new ArrayList<>();
            remainingPos.add(newPosIndex);
            for (int posNum = 2; posNum < posToMerge.size(); posNum++) {
                //posList[posNum]-1 to fix index after merging first two positions
                remainingPos.add(posToMerge.get(posNum) - 1);
            }
            mergePositionContinuous(toMergeNow);
            mergePositionContinuous(remainingPos);
        }
    }

    //This function will merge the positions and then merge the energy matrix
    ///efficiently based on what positions are being merged
    //posToMerge contains all the positions to merges
    public void mergePositionRigid(ArrayList<Integer> posToMerge) {
        //newPosList contains all the info about our new merged position
        //In this example newPosList will be [[0],[1,3],[2],[4]]
        ArrayList<ArrayList<Integer>> newPosList = getNewPosList(posToMerge);
        //get the index in the list that has multiple positions
        int newPosIndex = getNewPosIndex(newPosList);

        //If the number of positions to merge is 2, we can merge
        if (posToMerge.size() == 2) {
            //Create temp prune Mat to decide which RCs to merge
            PruningMatrix newPruneMat = mergePruneMat(newPosList);
            confSpaceSuper.mergePosition(newPosList, newPruneMat);
            this.mergeRigidEnergyMatrix(newPosList, newPruneMat);

        } //Otherwise, we recurse by merging the first two positoin
        else {
            ArrayList<Integer> toMergeNow = new ArrayList<>();
            toMergeNow.add(posToMerge.get(0));
            toMergeNow.add(posToMerge.get(1));
            ArrayList<Integer> remainingPos = new ArrayList<>();
            remainingPos.add(newPosIndex);
            for (int posNum = 2; posNum < posToMerge.size(); posNum++) {
                //posList[posNum]-1 to fix index after merging first two positions
                remainingPos.add(posToMerge.get(posNum) - 1);
            }
            mergePositionRigid(toMergeNow);
            mergePositionRigid(remainingPos);
        }
    }

    private PruningMatrix mergePruneMat(ArrayList<ArrayList<Integer>> newPosList) {
        int[] numRCPerNewPos = this.getNumRCPerNewPos(newPosList);
        PruningMatrix newPruneMat = new PruningMatrix(numRCPerNewPos.length, numRCPerNewPos, this.pruneMat.getPruningInterval());

        ArrayList<ArrayList<Boolean>> oneBody = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<ArrayList<Boolean>>>> twoBody = new ArrayList<>();

        for (int newPosIndex1 = 0; newPosIndex1 < newPosList.size(); newPosIndex1++) {
            ArrayList<Integer> newPos1 = newPosList.get(newPosIndex1);
            ArrayList<Boolean> oneBodyPrune = getMergedOneBodyPrune(newPos1);
            oneBody.add(oneBodyPrune);

            ArrayList<ArrayList<ArrayList<Boolean>>> pairwiseAtPos = new ArrayList<>();
            for (int newPosIndex2 = 0; newPosIndex2 < newPosIndex1; newPosIndex2++) {
                ArrayList<Integer> newPos2 = newPosList.get(newPosIndex2);
                ArrayList<ArrayList<Boolean>> pairwiseAtPair = getMergedTwoBodyPrune(newPos1, newPos2);
                pairwiseAtPos.add(pairwiseAtPair);
            }
            pairwiseAtPos.trimToSize();
            twoBody.add(pairwiseAtPos);
        }
        oneBody.trimToSize();
        twoBody.trimToSize();
        newPruneMat.oneBody = oneBody;
        newPruneMat.pairwise = twoBody;
        return newPruneMat;
    }

    private int[] getNumRCPerNewPos(ArrayList<ArrayList<Integer>> newPosList) {
        int[] originalNumRCPerPos = this.confSpaceSuper.getNumRCsAtPos();
        int[] numRCPerNewPos = new int[newPosList.size()];
        for (int newPosNum = 0; newPosNum < newPosList.size(); newPosNum++) {
            ArrayList<Integer> newPos = newPosList.get(newPosNum);
            if (newPos.size() == 1) {
                int oldPosNum = newPos.get(0);
                numRCPerNewPos[newPosNum] = originalNumRCPerPos[oldPosNum];
            } else {//newPos.size()==2
                int oldPosNum1 = newPos.get(0);
                int oldPosNum2 = newPos.get(1);
                numRCPerNewPos[newPosNum] = originalNumRCPerPos[oldPosNum1] * originalNumRCPerPos[oldPosNum2];
            }
        }
        return numRCPerNewPos;
    }

    private ArrayList<Boolean> getMergedOneBodyPrune(ArrayList<Integer> newPosList) {
        ArrayList<Boolean> oneBodyPrune = new ArrayList<>();
        if (newPosList.size() == 1) {
            int originalPosNum = newPosList.get(0);
            oneBodyPrune = this.pruneMat.oneBody.get(originalPosNum);
        } else {//newPosList.size()==2
            //iterate over larger position first
            int pos1 = Math.max(newPosList.get(0), newPosList.get(1));
            int numRCPos1 = this.emat.numRCsAtPos(pos1);
            int pos2 = Math.min(newPosList.get(0), newPosList.get(1));
            int numRCPos2 = this.emat.numRCsAtPos(pos2);
            for (int rc1 = 0; rc1 < numRCPos1; rc1++) {
                for (int rc2 = 0; rc2 < numRCPos2; rc2++) {
                    SuperRCTuple tup = new SuperRCTuple(pos1, rc1, pos2, rc2);
                    oneBodyPrune.add(this.pruneMat.isPruned(tup));
                }
            }
        }
        return oneBodyPrune;
    }

    private ArrayList<ArrayList<Boolean>> getMergedTwoBodyPrune(ArrayList<Integer> newPosList1, ArrayList<Integer> newPosList2) {
        //newPos1 > newPos2
        ArrayList<ArrayList<Boolean>> twoBodyPrune = new ArrayList<>();

        //if both are unmerged
        if (newPosList1.size() == 1 && newPosList2.size() == 1) {
            int posNum1 = Math.max(newPosList1.get(0), newPosList2.get(0));
            int posNum2 = Math.min(newPosList1.get(0), newPosList2.get(0));
            twoBodyPrune = this.pruneMat.pairwise.get(posNum1).get(posNum2);
        } else if (newPosList1.size() == 2 && newPosList2.size() == 1) {
            int oldPos1_1 = Math.max(newPosList1.get(0), newPosList1.get(1));
            int numRC1_1 = this.emat.numRCsAtPos(oldPos1_1);
            int oldPos1_2 = Math.min(newPosList1.get(0), newPosList1.get(1));
            int numRC1_2 = this.emat.numRCsAtPos(oldPos1_2);

            int oldPos2 = newPosList2.get(0);
            int numRC2 = this.emat.numRCsAtPos(oldPos2);

            for (int rc1_1 = 0; rc1_1 < numRC1_1; rc1_1++) {
                for (int rc1_2 = 0; rc1_2 < numRC1_2; rc1_2++) {
                    ArrayList<Boolean> pairwisePrune = new ArrayList<>();
                    for (int rc2 = 0; rc2 < numRC2; rc2++) {
                        ArrayList<Integer> posList = new ArrayList<>();
                        ArrayList<Integer> rcList = new ArrayList<>();
                        posList.add(oldPos1_1);
                        posList.add(oldPos1_2);
                        posList.add(oldPos2);
                        rcList.add(rc1_1);
                        rcList.add(rc1_2);
                        rcList.add(rc2);
                        SuperRCTuple tup = new SuperRCTuple(posList, rcList);
                        pairwisePrune.add(this.pruneMat.isPruned(tup));
                    }
                }
            }
        } else if (newPosList1.size() == 1 && newPosList2.size() == 2) {
            int oldPos1 = newPosList1.get(0);
            int numRC1 = this.emat.numRCsAtPos(oldPos1);

            int oldPos2_1 = Math.max(newPosList2.get(0), newPosList2.get(1));
            int numRC2_1 = this.emat.numRCsAtPos(oldPos2_1);
            int oldPos2_2 = Math.min(newPosList2.get(0), newPosList2.get(1));
            int numRC2_2 = this.emat.numRCsAtPos(oldPos2_2);

            for (int rc1 = 0; rc1 < numRC1; rc1++) {
                ArrayList<Boolean> pairwisePrune = new ArrayList<>();
                for (int rc2_1 = 0; rc2_1 < numRC2_1; rc2_1++) {
                    for (int rc2_2 = 0; rc2_2 < numRC2_2; rc2_2++) {
                        ArrayList<Integer> posList = new ArrayList<>();
                        ArrayList<Integer> rcList = new ArrayList<>();
                        posList.add(oldPos1);
                        posList.add(oldPos2_1);
                        posList.add(oldPos2_2);
                        rcList.add(rc1);
                        rcList.add(rc2_1);
                        rcList.add(rc2_2);
                        SuperRCTuple tup = new SuperRCTuple(posList, rcList);
                        pairwisePrune.add(this.pruneMat.isPruned(tup));
                    }
                }
            }
        } else {
            throw new RuntimeException("Cannot Merge twoBodyPrune: CAN ONLY MERGE TWO POSITIONS IN RIGID CASE!");
        }
        twoBodyPrune.trimToSize();
        return twoBodyPrune;
    }

    private ArrayList<ArrayList<Integer>> getNewPosList(ArrayList<Integer> posToMerge) {
        //newPosList contains all the info about our new merged position
        //In this example newPosList will be [[0],[1,3],[2],[4]]
        ArrayList<ArrayList<Integer>> newPosList = new ArrayList<>();

        //index into newPostList of new merged position
        //if newPosList is [[0],[1,3],[2],[4]], then newPosIndex = 1;
        int newPosIndex = -1;
        //new merge position consisting of multiple old positions
        ArrayList<Integer> newPos = new ArrayList<>();

        //iterate overall positions
        for (int i = 0; i < confSpaceSuper.numPos; i++) {
            //if this position is not going to be merged create a arraylist
            //that contains only this element
            if (!posToMerge.contains(i)) {
                ArrayList<Integer> pos = new ArrayList<>();
                pos.add(i);
                newPosList.add(pos);
            } //if this position is going to be merged
            else {
                //check if newPosIndex has been set
                if (newPosIndex > -1) {
                    //if it has been, then we add this position to the arraylist
                    //at the index of newPosList specified by newPosIndex
                    newPosList.get(newPosIndex).add(i);
                } else {
                    //if it has not been set then we set it and create the arraylist
                    //in this example i = 1 will reach this step.
                    newPosIndex = i;
                    ArrayList<Integer> mergePos = new ArrayList<>();
                    mergePos.add(i);
                    newPosList.add(mergePos);
                }
            }
        }
        return newPosList;
    }

    private int getNewPosIndex(ArrayList<ArrayList<Integer>> newPosList) {
        int index = -1;
        for (int i = 0; i < newPosList.size(); i++) {
            if (newPosList.get(i).size() > 1) {
                index = i;
            }
        }
        if (index == -1) {
            throw new RuntimeException("ERROR: There must be more than one positions to merge");
        }
        return index;
    }

    private void mergeContinuousEnergyMatrix(ArrayList<ArrayList<Integer>> newPosList, PruningMatrix pruneMat) {
        //intialize newEmat using info in newPosListCopy
        int[] numRCPerNewPos = this.confSpaceSuper.getNumRCsAtPos();
        EnergyMatrix newEmat = new EnergyMatrix(newPosList.size(), numRCPerNewPos, Double.POSITIVE_INFINITY);
        //one-body
        ArrayList<ArrayList<Double>> oneBody = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> twoBody = new ArrayList<>();

        for (int newPosIndex1 = 0; newPosIndex1 < newPosList.size(); newPosIndex1++) {
            ArrayList<Integer> newPos1 = newPosList.get(newPosIndex1);
            ArrayList<Double> oneBodyE = getMergedOneBodyContinuousEnergy(newPos1, newPosIndex1);
            oneBody.add(oneBodyE);

            ArrayList<ArrayList<ArrayList<Double>>> pairwiseAtPos = new ArrayList<>();
            for (int newPosIndex2 = 0; newPosIndex2 < newPosIndex1; newPosIndex2++) {
                ArrayList<Integer> newPos2 = newPosList.get(newPosIndex2);
                ArrayList<ArrayList<Double>> pairwiseAtPair = getMergedTwoBodyContinuousEnergy(newPos1, newPosIndex1, newPos2, newPosIndex2);
                pairwiseAtPos.add(pairwiseAtPair);
            }
            pairwiseAtPos.trimToSize();
            twoBody.add(pairwiseAtPos);
        }
        oneBody.trimToSize();
        twoBody.trimToSize();
        newEmat.oneBody = oneBody;
        newEmat.pairwise = twoBody;
        this.emat = newEmat;
        ObjectIO.writeObject(this.emat, name + ".EMAT_MERGED.dat");

    }

    private ArrayList<Double> getMergedOneBodyContinuousEnergy(ArrayList<Integer> newPos, int newPosNum) {
        ArrayList<Double> oneBodyE;
        //If this position is not being merged
        if (newPos.size() == 1) {
            oneBodyE = this.emat.oneBody.get(newPos.get(0));
        } else {
            TermECalculatorSuper termE = new TermECalculatorSuper(confSpaceSuper, shellResidues, useEPIC, useEPIC, null, epicSettings, false, newPosNum);
            Object oneBodyTerm = termE.doCalculation();
            oneBodyE = (ArrayList<Double>) oneBodyTerm;
        }
        return oneBodyE;
    }

    private ArrayList<ArrayList<Double>> getMergedTwoBodyContinuousEnergy(ArrayList<Integer> posNumList1, int newPosNum1, ArrayList<Integer> posNumList2, int newPosNum2) {
        ArrayList<ArrayList<Double>> twoBodyE;
        //
        //if either of the positions are merged
        if (posNumList1.size() > 1 || posNumList2.size() > 1) {
            TermECalculatorSuper termE = new TermECalculatorSuper(confSpaceSuper, shellResidues, useEPIC, useEPIC, null, epicSettings, false, newPosNum1, newPosNum2);
            Object twoBodyTerm = termE.doCalculation();
            twoBodyE = (ArrayList<ArrayList<Double>>) twoBodyTerm;
        } else {
            int originalPosNum1 = posNumList1.get(0);
            int originalPosNum2 = posNumList2.get(0);
            if (originalPosNum1 > originalPosNum2) {
                twoBodyE = this.emat.pairwise.get(originalPosNum1).get(originalPosNum2);
            } else {
                twoBodyE = this.emat.pairwise.get(originalPosNum2).get(originalPosNum1);
            }
        }
        return twoBodyE;
    }

    private void mergeRigidEnergyMatrix(ArrayList<ArrayList<Integer>> newPosList, PruningMatrix newPruneMat) {
        //intialize newEmat using info in newPosListCopy
        int[] numRCPerNewPos = this.confSpaceSuper.getNumRCsAtPos();
        EnergyMatrix newEmat = new EnergyMatrix(newPosList.size(), numRCPerNewPos, Double.POSITIVE_INFINITY);
        newEmat.setConstTerm(this.emat.getConstTerm());
        //one-body
        ArrayList<ArrayList<Double>> oneBody = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> twoBody = new ArrayList<>();

        for (int newPosIndex1 = 0; newPosIndex1 < newPosList.size(); newPosIndex1++) {
            ArrayList<Integer> newPos1 = newPosList.get(newPosIndex1);
            ArrayList<Double> oneBodyE = getMergedOneBodyRigidEnergy(newPos1, newPosIndex1, newPruneMat);
            oneBody.add(oneBodyE);

            ArrayList<ArrayList<ArrayList<Double>>> pairwiseAtPos = new ArrayList<>();
            for (int newPosIndex2 = 0; newPosIndex2 < newPosIndex1; newPosIndex2++) {
                ArrayList<Integer> newPos2 = newPosList.get(newPosIndex2);
                ArrayList<ArrayList<Double>> pairwiseAtPair = getMergedTwoBodyRigidEnergy(newPos1, newPosIndex1, newPos2, newPosIndex2, newPruneMat);
                pairwiseAtPos.add(pairwiseAtPair);
            }
            pairwiseAtPos.trimToSize();
            twoBody.add(pairwiseAtPos);
        }
        oneBody.trimToSize();
        twoBody.trimToSize();
        newEmat.oneBody = oneBody;
        newEmat.pairwise = twoBody;
        this.emat = newEmat;
        ObjectIO.writeObject(this.emat, name + ".EMAT_MERGED.dat");
    }

    private ArrayList<Double> getMergedOneBodyRigidEnergy(ArrayList<Integer> posNumList, int newPosNum, PruningMatrix pruneMat) {
        ArrayList<Double> oneBodyE = new ArrayList<>();
        if (posNumList.size() == 1) {
            int posNum = posNumList.get(0);
            oneBodyE = this.emat.oneBody.get(posNum);
        } else if (posNumList.size() == 2) {
            int pos1 = Math.max(posNumList.get(0), posNumList.get(1));
            int numRCPos1 = this.emat.numRCsAtPos(pos1);
            int pos2 = Math.min(posNumList.get(0), posNumList.get(1));
            int numRCPos2 = this.emat.numRCsAtPos(pos2);

            int totalRCNum = 0;
            //energy = oneBody(RC1) + oneBody(RC2) + pairwise(RC1,RC2)
            for (int rc1 = 0; rc1 < numRCPos1; rc1++) {
                double oldOneBodyE_rc1 = this.emat.getOneBody(pos1, rc1);
                for (int rc2 = 0; rc2 < numRCPos2; rc2++) {
                    double oldONeBodyE_rc2 = this.emat.getOneBody(pos2, rc2);
                    double oldPairwiseE = this.emat.getPairwise(pos1, rc1, pos2, rc2);
                    double newE = oldOneBodyE_rc1 + oldONeBodyE_rc2 + oldPairwiseE;
                    SuperRCTuple tup = new SuperRCTuple(newPosNum, totalRCNum);
                    //only add rotamer if it is not pruned (or parametrically incompatible)
                    if (!pruneMat.isPruned(tup)) {
                        oneBodyE.add(newE);
                    }
                    totalRCNum++;
                }
            }
        } else {
            throw new RuntimeException("Can't merge oneBodyE: CAN ONLY MERGE TWO POSITIONS IN RIGID CASE!");
        }

        oneBodyE.trimToSize();
        return oneBodyE;
    }

    private ArrayList<ArrayList<Double>> getMergedTwoBodyRigidEnergy(ArrayList<Integer> posNumList1, int newPosNum1, ArrayList<Integer> posNumList2, int newPosNum2, PruningMatrix pruneMat) {
        ArrayList<ArrayList<Double>> twoBodyE = new ArrayList<>();
        //posNum1 > posNum 2

        //if both are unMerged positions, then we use the old pairwiseE
        if (posNumList1.size() == 1 && posNumList2.size() == 1) {
            int posNum1 = posNumList1.get(0);
            int posNum2 = posNumList2.get(0);
            //posNum2 should be greater than posNum1
            twoBodyE = this.emat.pairwise.get(posNum1).get(posNum2);
        } else if (posNumList1.size() == 2 && posNumList2.size() == 1) {
            int newPos1_1 = posNumList1.get(0);
            int numRC1_1 = this.emat.numRCsAtPos(newPos1_1);
            int newPos1_2 = posNumList1.get(1);
            int numRC1_2 = this.emat.numRCsAtPos(newPos1_2);

            int posNum2 = posNumList2.get(0);
            int numRC2 = this.emat.numRCsAtPos(posNum2);

            int totalNewRCNum = 0;
            for (int rc1_1 = 0; rc1_1 < numRC1_1; rc1_1++) {
                for (int rc1_2 = 0; rc1_2 < numRC1_2; rc1_2++) {
                    SuperRCTuple tup = new SuperRCTuple(newPosNum1, totalNewRCNum);
                    //if this new SuperRC is pruned, it will not be added
                    if (!pruneMat.isPruned(tup)) {
                        ArrayList<Double> pairwiseE = new ArrayList<>();
                        for (int rc2 = 0; rc2 < numRC2; rc2++) {
                            double pairwise1_2 = this.emat.getPairwise(newPos1_1, rc1_1, posNum2, rc2);
                            double pairwise2_2 = this.emat.getPairwise(newPos1_2, rc1_2, posNum2, rc2);
                            double newPairwiseE = pairwise1_2 + pairwise2_2;
                            pairwiseE.add(newPairwiseE);
                        }
                        pairwiseE.trimToSize();
                        twoBodyE.add(pairwiseE);
                    }
                    totalNewRCNum++;
                }
            }
        } else if (posNumList1.size()
                == 1 && posNumList2.size() == 2) {
            int newPos2_1 = posNumList2.get(0);
            int numRC2_1 = this.emat.numRCsAtPos(newPos2_1);
            int newPos2_2 = posNumList2.get(1);
            int numRC2_2 = this.emat.numRCsAtPos(newPos2_2);

            int posNum1 = posNumList1.get(0);
            int numRC1 = this.emat.numRCsAtPos(posNum1);

            for (int rc1 = 0; rc1 < numRC1; rc1++) {
                ArrayList<Double> pairwiseE = new ArrayList<>();
                int totalNewRCNum = 0;
                for (int rc2_1 = 0; rc2_1 < numRC2_1; rc2_1++) {
                    for (int rc2_2 = 0; rc2_2 < numRC2_2; rc2_2++) {
                        SuperRCTuple tup = new SuperRCTuple(newPosNum2, totalNewRCNum);
                        //if new SuperRC is pruned, we will not add it
                        if (!pruneMat.isPruned(tup)) {
                            double pairwise1_2 = this.emat.getPairwise(posNum1, rc1, newPos2_1, rc2_1);
                            double pairwise2_2 = this.emat.getPairwise(posNum1, rc1, newPos2_2, rc2_2);
                            double newPairwiseE = pairwise1_2 + pairwise2_2;
                            pairwiseE.add(newPairwiseE);
                        }
                        totalNewRCNum++;
                    }
                }
                pairwiseE.trimToSize();
                twoBodyE.add(pairwiseE);
            }
        } else {
            throw new RuntimeException("Cannot Merge twoBodyE: CAN ONLY MERGE TWO POSITIONS IN RIGID CASE!");
        }
        twoBodyE.trimToSize();
        return twoBodyE;
    }

//load the specified matrix; if the right file isn't available then compute and store it
    private void loadMatrix(MatrixType type) {

        String matrixFileName = name + "." + type.name() + ".dat";
        TupleMatrix matrix = new TupleMatrix();
        //matrix file names are determined by the name of the search problem
        if (!loadMatrixFromFile(type, matrixFileName)) {
            matrix = calcMatrix(type);
            ObjectIO.writeObject(matrix, matrixFileName);
            loadMatrixFromFile(type, matrixFileName);
        }

    }

    public void loadMergedEnergyMatrix() {
        MatrixType type = MatrixType.EMAT;
        String matrixFileName = name + ".EMAT_MERGED.dat";
        loadMatrixFromFile(type, matrixFileName);
    }

    //compute the matrix of the specified type
    private TupleMatrix calcMatrix(MatrixType type) {

        if (type == MatrixType.EMAT) {
            EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpaceSuper, shellResidues);
            emCalc.calcPEM();
            return emCalc.getEMatrix();
        } else if (type == MatrixType.EPICMAT) {
            throw new RuntimeException("TupExpansion and EPIC not currently supported with SearchProblemSuper and super-RCs");
            /*
             EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpaceSuper, shellResidues,
             pruneMat, epicSettings);
             emCalc.calcPEM();
             return emCalc.getEPICMatrix();
             */
        } else {
            throw new RuntimeException("ONLY MatrixType.EMAT is currently supported in SearchProblemSuper and super-RCs");

            /*
             //need to calculate a tuple-expansion matrix

             double errorThresh = 0.01;

             ConfETupleExpander expander = new ConfETupleExpander(this);//make a tuple expander
             TupleEnumerator tupEnum = new TupleEnumerator(pruneMat, emat, confSpace.numPos);
             TupExpChooser chooser = new TupExpChooser(expander, tupEnum);//make a chooser to choose what tuples will be in the expansion

             double curResid = chooser.calcPairwiseExpansion();//start simple...

             if (curResid > errorThresh) {//go to triples if needed
             System.out.println("EXPANDING PAIRWISE EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (2 PARTNERS)...");
             curResid = chooser.calcExpansionResTriples(2);
             }
             if (curResid > errorThresh) {//go to 5 partners if still need better resid...
             System.out.println("EXPANDING EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (5 PARTNERS)...");
             curResid = chooser.calcExpansionResTriples(5);
             }
             if (curResid > errorThresh) {
             System.out.println("WARNING: Desired LUTE residual threshold "
             + errorThresh + " not reached; best=" + curResid);
             }

             return expander.getEnergyMatrix();//get the final energy matrix from the chosen expansion
             }
             */
        }
    }

    boolean loadMatrixFromFile(MatrixType type, String matrixFileName) {
        //try loading the specified matrix from a file
        //return true if successful, false if not, in which case we'll have to compute it
        //also if the matrix's pruning interval is too low, it may be missing some RCs
        //that are unpruned at our current pruningInterval, so we have to recompute
        Object matrixFromFile = ObjectIO.readObject(matrixFileName, true);

        if (type == MatrixType.EMAT) {
            emat = (EnergyMatrix) matrixFromFile;
        } else if (type == MatrixType.EPICMAT) {
            epicMat = (EPICMatrix) matrixFromFile;
        } else //tup-exp
        {
            tupExpEMat = (EnergyMatrix) matrixFromFile;
        }

        if (matrixFromFile == null)//unsuccessful loading leaves null emat
        {
            return false;
        }

        //check pruning interval.  Current interval is in pruneMat if we have pruned already;
        //if not then we need a matrix with infinite pruning interval (valid for all RCs).
        double matrixPruningInterval = ((TupleMatrix) matrixFromFile).getPruningInterval();

        if (matrixPruningInterval == Double.POSITIVE_INFINITY)//definitely valid
        {
            return true;
        } else {
            //excludes some RCs...check against pruneMat pruning interval
            if (pruneMat == null) {
                throw new RuntimeException("ERROR: Trying to load pruning-dependent tuple matrix"
                        + "(EPIC or tup-exp) but haven't pruned yet");
            }

            return (matrixPruningInterval >= pruneMat.getPruningInterval());
        }
    }
}
