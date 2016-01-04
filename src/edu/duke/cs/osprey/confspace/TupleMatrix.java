/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.tools.ObjectIO;
import java.io.Serializable;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;
import java.util.stream.Collectors;

/**
 *
 * @author mhall44
 */
public class TupleMatrix<T> implements Serializable {
    //We will need "matrices" of quantities defined
    //for example, the energy matrix (T=Double) stores single, pairwise, and higher-order energies
    //and we'll also have pruning (T=Boolean) and EPIC (T=EPoly) matrices
    //we'll store things as ArrayLists to make it easier to merge residues, partition RCs, etc.
    //and also to facilitate generics

    //note: tuples are sets not ordered pairs, i.e. E(i_r,j_s) = E(j_s,i_r), and pruning (i_r,j_s) means pruning (j_s,i_r)
    public ArrayList<ArrayList<ArrayList<ArrayList<T>>>> pairwise;//pairwise energies.  Set up 4D to save space
    //indices: res1, res2, RC1, RC2 where res1>res2
    public ArrayList<ArrayList<T>> oneBody;//intra+shell

    public ArrayList<ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>>> higherTerms;//look up higher terms by pair
    //same indices as pairwise
    //can be null if no interactions

    //maybe separate intra too?
    //The above all use RCs that are indexed by residues
    //the following lists indicate exactly what RCs those are
    //private because changing it may require re-indexing everything else
    //METHODS TO ADD RCs?
    //private ArrayList<ArrayList<RC>> RCList;
    /*
     //reverse lookup: for each residue, a lookup first by AA type, then by 
     //Return the residue-based RC number
     ArrayList<TreeMap<String,TreeMap<Integer,
     */
    double pruningInterval;//This matrix needs to hold entries for all RCs
    //that cannot be pruned with the specified pruning interval (Ew + Ival)
    //i.e. the matrix must describe all conformations within pruningInterval 
    //of the lowest pairwise lower bound

    T defaultHigherInteraction;//We only mark sparse higher interactions;
    //if unmarked we assume this value (e.g., 0 for energy, false for pruning)

    public TupleMatrix(T defaultHigherInteraction) {
        //no allocation (for overriding)
        this.defaultHigherInteraction = defaultHigherInteraction;
    }

    public TupleMatrix() {
    }

    public TupleMatrix(ConfSpace cSpace, double pruningInterval, T defaultHigherInteraction) {
        //allocate the matrix based on the provided conformational space
        init(cSpace.numPos, cSpace.getNumRCsAtPos(), pruningInterval, defaultHigherInteraction);
    }

    //HMN: Constructuro for input = ConfSpaceSuper
    ///This is unnecessary if ConfSpaceSuper inherits from ConfSpace, but we may
    ///remove this inheritance later
    public TupleMatrix(ConfSpaceSuper cSpace, double pruningInterval, T defaultHigherInteraction) {
        init(cSpace.numPos, cSpace.getNumRCsAtPos(), pruningInterval, defaultHigherInteraction);
    }

    public TupleMatrix(int numPos, int[] numAllowedAtPos, double pruningInterval, T defaultHigherInteraction) {
        //allocate the matrix based on the provided conformational space size
        //also specify what pruningInterval it's valid up to
        init(numPos, numAllowedAtPos, pruningInterval, defaultHigherInteraction);
    }

    public TupleMatrix(TupleMatrix<T> tupMat) {
        this.oneBody = tupMat.oneBody;
        this.pairwise = tupMat.pairwise;
        this.higherTerms = tupMat.higherTerms;
        this.pruningInterval = tupMat.pruningInterval;
        this.defaultHigherInteraction = tupMat.defaultHigherInteraction;
    }

    private void init(int numPos, int[] numAllowedAtPos, double pruningInterval, T defaultHigherInteraction) {

        this.pruningInterval = pruningInterval;
        this.defaultHigherInteraction = defaultHigherInteraction;

        oneBody = new ArrayList<>();
        pairwise = new ArrayList<>();

        higherTerms = new ArrayList<>();//preallocate these too, but all null for now (no higher-order terms yet)

        for (int pos = 0; pos < numPos; pos++) {

            int numRCs = numAllowedAtPos[pos];

            //preallocate oneBody for this position
            ArrayList<T> oneBodyAtPos = new ArrayList<>();
            for (int rc = 0; rc < numRCs; rc++)//preallocate oneBody
            {
                oneBodyAtPos.add(null);
            }

            oneBodyAtPos.trimToSize();//we may need to save space so we'll trim everything to size
            oneBody.add(oneBodyAtPos);

            ArrayList<ArrayList<ArrayList<T>>> pairwiseAtPos = new ArrayList<>();
            //may want to leave some pairs of positions null if negligible interaction expected...
            //handle later though

            ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>> higherOrderAtPos = new ArrayList<>();

            for (int pos2 = 0; pos2 < pos; pos2++) {

                int numRCs2 = numAllowedAtPos[pos2];

                ArrayList<ArrayList<T>> pairwiseAtPair = new ArrayList<>();
                ArrayList<ArrayList<HigherTupleFinder<T>>> higherOrderAtPair = new ArrayList<>();

                for (int rc = 0; rc < numRCs; rc++) {
                    ArrayList<T> pairwiseAtRC = new ArrayList<>();
                    ArrayList<HigherTupleFinder<T>> higherOrderAtRC = new ArrayList<>();

                    for (int rc2 = 0; rc2 < numRCs2; rc2++) {
                        pairwiseAtRC.add(null);
                        higherOrderAtRC.add(null);
                    }

                    pairwiseAtRC.trimToSize();
                    pairwiseAtPair.add(pairwiseAtRC);

                    higherOrderAtRC.trimToSize();
                    higherOrderAtPair.add(higherOrderAtRC);
                }

                pairwiseAtPair.trimToSize();
                pairwiseAtPos.add(pairwiseAtPair);

                higherOrderAtPair.trimToSize();
                higherOrderAtPos.add(higherOrderAtPair);
            }

            pairwiseAtPos.trimToSize();
            pairwise.add(pairwiseAtPos);

            higherOrderAtPos.trimToSize();
            higherTerms.add(higherOrderAtPos);
        }

        oneBody.trimToSize();
        pairwise.trimToSize();
        higherTerms.trimToSize();
    }

    public T getPairwise(int res1, int index1, int res2, int index2) {
        //working with residue-specific RC indices directly.  
        if (res1 > res2) {
            return pairwise.get(res1).get(res2).get(index1).get(index2);
        } else {
            return pairwise.get(res2).get(res1).get(index2).get(index1);
        }
    }

    public T getOneBody(int res, int index) {
        return oneBody.get(res).get(index);
    }

    public void setPairwise(int res1, int index1, int res2, int index2, T val) {

        if (res1 > res2) {
            pairwise.get(res1).get(res2).get(index1).set(index2, val);
        } else {
            pairwise.get(res2).get(res1).get(index2).set(index1, val);
        }
    }

    public void setOneBody(int res, int index, T val) {
        oneBody.get(res).set(index, val);
    }

    public int numRCsAtPos(int pos) {
        return oneBody.get(pos).size();
    }

    public int numPos() {
        return oneBody.size();
    }

    public void setTupleValue(RCTuple tup, T val) {
        //assign the given value to the specified RC tuple
        int tupSize = tup.pos.size();

        if (tupSize == 1)//just a one-body quantity
        {
            setOneBody(tup.pos.get(0), tup.RCs.get(0), val);
        } else if (tupSize == 2)//two-body
        {
            setPairwise(tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1), val);
        } else if (tupSize > 2) {//higher-order
            setHigherOrder(tup, val);
        } else {
            throw new UnsupportedOperationException("ERROR: Not supporting tuple size " + tupSize);
        }
    }

    public void setTupleValue(SuperRCTuple tup, T val) {
        //assign the given value to the specified RC tuple
        int tupSize = tup.pos.size();
        if (tupSize == 1)//just a one-body quantity
        {
            setOneBody(tup.pos.get(0), tup.superRCs.get(0), val);
        } else if (tupSize == 2)//two-body
        {
            setPairwise(tup.pos.get(0), tup.superRCs.get(0), tup.pos.get(1), tup.superRCs.get(1), val);
        } else if (tupSize > 2) {//higher-order
            setHigherOrder(tup, val);
        } else {
            throw new UnsupportedOperationException("ERROR: Not supporting tuple size " + tupSize);
        }
    }

    public void setHigherOrder(RCTuple tup, T val) {
        //set a higher-order term
        //we need all pairs contained in tup to know about it

        //loop over pairs
        for (int index1 = 0; index1 < tup.pos.size(); index1++) {
            for (int index2 = 0; index2 < index1; index2++) {

                int pos1 = tup.pos.get(index1);
                int rc1 = tup.RCs.get(index1);
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);

                //put tup into the HigherTupleFinder for this pair
                HigherTupleFinder<T> htf = getHigherOrderTerms(pos1, rc1, pos2, rc2);

                //create a HigherTupleFinder if there is none yet
                if (htf == null) {
                    htf = new HigherTupleFinder(defaultHigherInteraction);

                    if (pos1 > pos2) {
                        higherTerms.get(pos1).get(pos2).get(rc1).set(rc2, htf);
                    } else {
                        higherTerms.get(pos2).get(pos1).get(rc2).set(rc1, htf);
                    }
                }

                RCTuple subTup = tup.subtractMember(index1).subtractMember(index2);

                htf.setInteraction(subTup, val);
            }
        }

    }

    public void setHigherOrder(SuperRCTuple tup, T val) {
        //set a higher-order term
        //we need all pairs contained in tup to know about it

        //loop over pairs
        for (int index1 = 0; index1 < tup.pos.size(); index1++) {
            for (int index2 = 0; index2 < index1; index2++) {

                int pos1 = tup.pos.get(index1);
                int rc1 = tup.superRCs.get(index1);
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.superRCs.get(index2);

                //put tup into the HigherTupleFinder for this pair
                HigherTupleFinder<T> htf = getHigherOrderTerms(pos1, rc1, pos2, rc2);

                //create a HigherTupleFinder if there is none yet
                if (htf == null) {
                    htf = new HigherTupleFinder(defaultHigherInteraction);

                    if (pos1 > pos2) {
                        higherTerms.get(pos1).get(pos2).get(rc1).set(rc2, htf);
                    } else {
                        higherTerms.get(pos2).get(pos1).get(rc2).set(rc1, htf);
                    }
                }

                SuperRCTuple subTup = tup.subtractMember(index1).subtractMember(index2);

                htf.setInteraction(subTup, val);
            }
        }

    }

    public double getPruningInterval() {
        return pruningInterval;
    }

    public HigherTupleFinder<T> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
        //working with residue-specific RC indices directly.  
        if (res1 > res2) {
            return higherTerms.get(res1).get(res2).get(index1).get(index2);
        } else {
            return higherTerms.get(res2).get(res1).get(index2).get(index1);
        }
    }

    //HMN: Useful for Partial Search Spaces
    public TupleMatrix<T> getSubsetMatrix(ArrayList<Integer> subsetOfPos) {
        Collections.sort(subsetOfPos);

        int numSubsetPos = subsetOfPos.size();
        int[] numRCsAtPos = new int[numSubsetPos];
        for (int pos = 0; pos < numSubsetPos; pos++) {
            int originalPosNum = subsetOfPos.get(pos);
            numRCsAtPos[pos] = this.oneBody.get(originalPosNum).size();
        }
        TupleMatrix<T> subsetTupMat = new TupleMatrix<>(numSubsetPos, numRCsAtPos, pruningInterval, defaultHigherInteraction);

        ArrayList<ArrayList<T>> newOneBody = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<ArrayList<T>>>> newTwoBody = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>>> newHigherOrder = new ArrayList<>();

        for (int i = 0; i < numSubsetPos; i++) {
            int originalPosNumI = subsetOfPos.get(i);
            ArrayList<T> oneBodyAtPosI = (ArrayList<T>) ObjectIO.deepCopy(this.oneBody.get(originalPosNumI));
            ArrayList<ArrayList<ArrayList<T>>> twoBodyAtPosI = new ArrayList<>();
            ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>> higherOrderAtPosI = new ArrayList<>();
            for (int j = 0; j < i; j++) {
                int originalPosNumJ = subsetOfPos.get(j);
                ArrayList<ArrayList<T>> twoBodyAtPosIPosJ = (ArrayList<ArrayList<T>>) ObjectIO.deepCopy(this.pairwise.get(originalPosNumI).get(originalPosNumJ));
                twoBodyAtPosI.add(twoBodyAtPosIPosJ);
                //Higher-Order Tuples
                ArrayList<ArrayList<HigherTupleFinder<T>>> higherOrderAtPosIPosJ = new ArrayList<>();
                for (int rcI = 0; rcI < numRCsAtPos[i]; rcI++) {
                    ArrayList<HigherTupleFinder<T>> higherOrderAtPosIPosJRCI = new ArrayList<>();
                    for (int rcJ = 0; rcJ < numRCsAtPos[j]; rcJ++) {
                        HigherTupleFinder<T> htf = getHigherOrderTerms(originalPosNumI, rcI, originalPosNumJ, rcJ);
                        if (htf != null) {
                            HigherTupleFinder<T> new_htf = htf.getHigherOrderTuplesWithinSubsetPos(subsetOfPos, numRCsAtPos);
                            higherOrderAtPosIPosJRCI.add(new_htf);

                        } else {
                            higherOrderAtPosIPosJRCI.add(null);
                        }
                    }
                    higherOrderAtPosIPosJ.add(higherOrderAtPosIPosJRCI);
                }
                higherOrderAtPosI.add(higherOrderAtPosIPosJ);
            }
            newOneBody.add(oneBodyAtPosI);
            newTwoBody.add(twoBodyAtPosI);
            newHigherOrder.add(higherOrderAtPosI);
        }
        subsetTupMat.oneBody = newOneBody;
        subsetTupMat.pairwise = newTwoBody;
        subsetTupMat.higherTerms = newHigherOrder;

        return subsetTupMat;
    }

    
    /*
    public void updateMatrixCrossTerms(TupleMatrix<T> tupMat, boolean[][] interactionGraph) {
        int numPos = this.oneBody.size();

        int[] numRCsPerPos = new int[numPos];
        for (int pos=0; pos<numPos; pos++){
            numRCsPerPos[pos] = this.oneBody.get(pos).size();
        }
        
        
        for (int i = 0; i < numPos; i++) {
            //Set all one-body terms to "zero"
            for (int rot=0; rot< numRCsPerPos[i]; rot++){
                this.setOneBody(i, rot, defaultHigherInteraction);
            }
            
            for (int j = 0; j < i; j++) {
                //If they don't interact, set all pairwise to "zero"
                //This is the defaultHigherInteraction (e.g. false, null, 0)
                if (!interactionGraph[i][j]) {
                    for (int rcI=0; rcI<numRCsPerPos[i]; rcI++){
                        for (int rcJ=0; rcJ<numRCsPerPos[j]; rcJ++){
                            //Set pairwise to "zero"
                            this.setPairwise(i, rcI, j, rcJ, defaultHigherInteraction);
                            //now we need to check if we should include higher-order terms
                            HigherTupleFinder<T> htf = getHigherOrderTerms(i, rcI, j, rcJ);
                            if (htf != null){
                                //Get all positions that interact with i or j
                                ArrayList<Integer> interactingPositions = new ArrayList<>();
                                for (int pos=0; pos<numPos; pos++){
                                    if (interactionGraph[i][pos] || interactionGraph[j][pos]){
                                        interactingPositions.add(pos);
                                    }
                                }
                                HigherTupleFinder<T> crossHtf = htf.getHigherOrderCrossTerms(interactingPositions, numRCsPerPos);
                                this.higherTerms.get(i).get(j).get(rcI).set(rcJ, crossHtf);
                            }
                        }
                    }
                }
            }
        }
    }
    */
}
