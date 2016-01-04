/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.tools.ObjectIO;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 *
 * This is an object associated to a particular tuple t (pair or higher) in a
 * TupleMatrix It lists any higher-tuples that include t and that have
 * non-trivial values ("trivial" meaning 0 for energies, unpruned for pruning)
 * Tuples of order size(t)+1 are held explicitly; we recursively list
 * HigherTupleFinder objects for higher orders if there are even bigger tuples
 *
 * @author mhall44
 */
public class HigherTupleFinder<T> implements Serializable {

    private ArrayList<Integer> interactingPos = new ArrayList<>();
    //flexible positions interacting with t in higher tuples

    //For sparseness purposes, the following lists are indexed by interacting position,
    //not by the usual position number
    //They map RC numbers (at the specified position) to the value of interest
    private ArrayList<TreeMap<Integer, T>> interactions = new ArrayList<>();
    private ArrayList<TreeMap<Integer, HigherTupleFinder<T>>> higher = new ArrayList<>();//next order up...

    //Note: the RC number could potentially be a code (e.g., -2) indicating a set of RCs,
    //not just one RC.
    //if so the TupleMatrix should store some indication of what these codes are
    T defaultInteraction;//If no interaction included explicitly, this is the default value
    //(e.g., 0 for energies, false for pruning)

    public HigherTupleFinder(T defaultInteraction) {
        //generally start empty
        this.defaultInteraction = defaultInteraction;
    }

    public ArrayList<Integer> getInteractingPos() {
        return interactingPos;
    }

    public T getInteraction(int posNum, int RCNum) {
        //get interaction of this tuple with the RC (posNum,RCNum)

        for (int interactionNum = 0; interactionNum < interactingPos.size(); interactionNum++) {
            //there should be only a few interacting positions so we can loop over them quickly

            if (interactingPos.get(interactionNum) == posNum) {

                if (interactions.get(interactionNum).containsKey(RCNum)) {
                    return interactions.get(interactionNum).get(RCNum);
                } else {
                    return defaultInteraction;
                }
            }
        }

        //if we get here no interaction (i.e., just default)
        return defaultInteraction;
    }

    public HigherTupleFinder<T> getHigherInteractions(int posNum, int RCNum) {
        //get higher-order interactions involving the super-tuple 
        //consisting of this tuple plus the RC (posNum,RCNum)

        for (int interactionNum = 0; interactionNum < interactingPos.size(); interactionNum++) {
            //there should be only a few interacting positions so we can loop over them quickly

            if (interactingPos.get(interactionNum) == posNum) {

                if (higher.get(interactionNum).containsKey(RCNum)) {
                    return higher.get(interactionNum).get(RCNum);
                } else {
                    return null;
                }
            }
        }

        //if we get here no higher interaction
        return null;
    }

    public void setInteraction(RCTuple tup, T val) {
        //set the interaction of this tuple with tup to the given value

        if (tup.pos.size() == 1) {//store interaction directly in this HigherTupleFinder
            int pos = tup.pos.get(0);
            int rc = tup.RCs.get(0);

            int posIndex = getPosIndex(pos);

            interactions.get(posIndex).put(rc, val);
        } else {//kick it up to the next level.  Make sure all sub-tuples of tup know about this interaction.  

            for (int index = 0; index < tup.pos.size(); index++) {
                int pos = tup.pos.get(index);
                int rc = tup.RCs.get(index);

                int posIndex = getPosIndex(pos);

                RCTuple subTup = tup.subtractMember(index);

                HigherTupleFinder<T> nextHTF = higher.get(posIndex).get(rc);

                if (nextHTF == null) {//allocate next-level HTF if not currently existent
                    nextHTF = new HigherTupleFinder<>(defaultInteraction);
                    higher.get(posIndex).put(rc, nextHTF);
                }

                nextHTF.setInteraction(subTup, val);
            }
        }
    }

    public void setInteraction(SuperRCTuple tup, T val) {
        //set the interaction of this tuple with tup to the given value

        if (tup.pos.size() == 1) {//store interaction directly in this HigherTupleFinder
            int pos = tup.pos.get(0);
            int superRC = tup.superRCs.get(0);

            int posIndex = getPosIndex(pos);

            interactions.get(posIndex).put(superRC, val);
        } else {//kick it up to the next level.  Make sure all sub-tuples of tup know about this interaction.  

            for (int index = 0; index < tup.pos.size(); index++) {
                int pos = tup.pos.get(index);
                int rc = tup.superRCs.get(index);

                int posIndex = getPosIndex(pos);

                SuperRCTuple subTup = tup.subtractMember(index);

                HigherTupleFinder<T> nextHTF = higher.get(posIndex).get(rc);

                if (nextHTF == null) {//allocate next-level HTF if not currently existent
                    nextHTF = new HigherTupleFinder<>(defaultInteraction);
                    higher.get(posIndex).put(rc, nextHTF);
                }

                nextHTF.setInteraction(subTup, val);
            }
        }
    }

    private int getPosIndex(int pos) {
        //find the index in interactingPos of pos
        //Create one if non-existent
        //This function is called when we need to set new interactions

        for (int interactionNum = 0; interactionNum < interactingPos.size(); interactionNum++) {
            //there should be only a few interacting positions so we can loop over them quickly

            if (interactingPos.get(interactionNum) == pos) {
                return interactionNum;
            }
        }

        //not found...create
        interactingPos.add(pos);
        interactions.add(new TreeMap<Integer, T>());
        higher.add(new TreeMap<Integer, HigherTupleFinder<T>>());
        return interactingPos.size() - 1;
    }

    public ArrayList<RCTuple> listInteractionsWithValue(T val) {
        //list tuples with the given interaction value
        //stored here

        ArrayList<RCTuple> ans = new ArrayList<>();

        //search recursively, recording tuples in decsending order of position
        recordInteractionsWithValue(val, ans, Integer.MAX_VALUE);

        return ans;
    }

    /**
     * Hunter: Given a subset of the total number of positions in a search
     * space, this method creates a higher-order tuple finder for all
     * higher-order interactions involving this tuple and the subset of
     * positions. This implementation is completely recursive
     *
     * @param subsetPos the subset of total positions
     * @return a new HigherTupleFinder containing higher-order tuples involving
     * this tuple and the subset of positions
     */
    public HigherTupleFinder<T> getHigherOrderTuplesWithinSubsetPos(ArrayList<Integer> subsetPos, int[] numRCsPerNewPos) {
        //Create a new HigherTupleFinder to be returned
        HigherTupleFinder<T> htf = new HigherTupleFinder<>(this.defaultInteraction);
        //For each interacting position
        for (int i = 0; i < this.interactingPos.size(); i++) {
            int originalPosNum = this.interactingPos.get(i);
            //is this position in the subset of positions we are interested in
            if (subsetPos.contains(originalPosNum)) {
                //our new search space has positions i=0,i=1,...i=subsetPos.size()
                //newPosNum is the index of the original posNum in subsetPos
                int newPosNum = subsetPos.indexOf(originalPosNum);
                ///iterate over the rc's for this position 
                for (int rc = 0; rc < numRCsPerNewPos[newPosNum]; rc++) {
                    //get the value of this interaction
                    RCTuple tup = new RCTuple(newPosNum, rc);
                    if (this.interactions.get(i).containsKey(rc)) {
                        T val = this.interactions.get(i).get(rc);

                        //add this interaction to our HigherTupleFinder
                        if (!htf.interactingPos.contains(newPosNum)) {
                            htf.interactingPos.add(newPosNum);
                            htf.interactions.add(new TreeMap<>());
                            htf.higher.add(new TreeMap<>());
                        }
                        int htfIndex = htf.interactingPos.indexOf(newPosNum);
                        htf.interactions.get(htfIndex).put(rc, val);
                    }
                    //Now we must check for higher-order
                    if (this.higher.get(i).containsKey(rc)) {
                        //In this case we have a higher-order interaction between 
                        //at least three positions in subsetPos
                        HigherTupleFinder<T> higherHTF = this.higher.get(i).get(rc);
                        //we recurse on the higher-order htf 
                        HigherTupleFinder<T> higherHTFsubset = higherHTF.getHigherOrderTuplesWithinSubsetPos(subsetPos, numRCsPerNewPos);
                        //If any higher-order interactions are found
                        if (higherHTFsubset != null) {
                            //... we add the higher-order tuple to htf
                            // first get into htf.higher
                            if (!htf.interactingPos.contains(newPosNum)) {
                                htf.interactingPos.add(newPosNum);
                                htf.interactions.add(new TreeMap<>());
                                htf.higher.add(new TreeMap<>());
                            }
                            int htfIndex = htf.interactingPos.indexOf(newPosNum);
                            htf.higher.get(htfIndex).put(rc, higherHTFsubset);
                        }
                    }
                }
            }
        }
        if (htf.interactingPos.isEmpty()) {
            return null;
        }
        return htf;
    }


    /**
     * This function gets higher-order terms for partial search spaces involving 
     * cross-term between subsets of residues
     * @param crossTermPos the set of positions that we want higher-order terms 
     * for this cross term
     * @param numRCsPerNewPos number of RC's per position
     * @return 
     */
    public HigherTupleFinder<T> getHigherOrderCrossTerms(ArrayList<Integer> crossTermPos, int[] numRCsPerNewPos) {
        HigherTupleFinder<T> htf = new HigherTupleFinder<>(this.defaultInteraction);
        for (int pos : this.interactingPos) {
            int currentPosIndex = this.interactingPos.indexOf(pos);

            //If this position is contained in our list of cross-terms, we add
            //all higher-order terms
            if (crossTermPos.contains(pos)) {
                //allocate space for this position in htf
                if (!htf.interactingPos.contains(pos)) {
                    htf.interactingPos.add(pos);
                    htf.interactions.add(new TreeMap<>());
                    htf.higher.add(new TreeMap<>());
                }
                //set interactions and higher
                int htfIndex = htf.interactingPos.indexOf(pos);
                this.interactions.set(htfIndex, (TreeMap<Integer, T>) ObjectIO.deepCopy(this.interactions.get(currentPosIndex)));
                this.higher.set(htfIndex, (TreeMap<Integer, HigherTupleFinder<T>>) ObjectIO.deepCopy(this.higher.get(currentPosIndex)));
            }
        }
        if (htf.interactingPos.isEmpty()){
            return null;
        }
        return htf;
    }

    public ArrayList<SuperRCTuple> listInteractionsWithValueSuper(T val) {
        //list tuples with the given interaction value
        //stored here

        ArrayList<SuperRCTuple> ans = new ArrayList<>();

        //search recursively, recording tuples in decsending order of position
        recordInteractionsWithValueSuper(val, ans, Integer.MAX_VALUE);

        return ans;
    }

    private void recordInteractionsWithValue(T val, ArrayList<RCTuple> tupList, int maxPos) {
        //Find tuples with all pos numbers less than maxPos that have interactions value val
        //add them to tupList

        for (int interactionNum = 0; interactionNum < interactingPos.size(); interactionNum++) {

            int pos = interactingPos.get(interactionNum);

            if (pos < maxPos) {

                for (int rc : interactions.get(interactionNum).keySet()) {

                    if (interactions.get(interactionNum).get(rc) == val) {
                        tupList.add(new RCTuple(pos, rc));
                    }
                }

                for (int rc : higher.get(interactionNum).keySet()) {

                    ArrayList<RCTuple> subTupList = new ArrayList<>();
                    higher.get(interactionNum).get(rc).recordInteractionsWithValue(val, subTupList, pos);
                    //the HigherTupleFinders in higher don't know about (pos,rc), so add that

                    for (RCTuple subTup : subTupList) {
                        subTup.pos.add(pos);
                        subTup.RCs.add(rc);
                        tupList.add(subTup);
                    }
                }
            }
        }
    }

    //Same as above recordInteractionsWithValue but with SuperRCTuple
    private void recordInteractionsWithValueSuper(T val, ArrayList<SuperRCTuple> tupList, int maxPos) {
        //Find tuples with all pos numbers less than maxPos that have interactions value val
        //add them to tupList

        for (int interactionNum = 0; interactionNum < interactingPos.size(); interactionNum++) {

            int pos = interactingPos.get(interactionNum);

            if (pos < maxPos) {

                for (int rc : interactions.get(interactionNum).keySet()) {

                    if (interactions.get(interactionNum).get(rc) == val) {
                        tupList.add(new SuperRCTuple(pos, rc));
                    }
                }

                for (int rc : higher.get(interactionNum).keySet()) {

                    ArrayList<SuperRCTuple> subTupList = new ArrayList<>();
                    higher.get(interactionNum).get(rc).recordInteractionsWithValueSuper(val, subTupList, pos);
                    //the HigherTupleFinders in higher don't know about (pos,rc), so add that

                    for (SuperRCTuple subTup : subTupList) {
                        subTup.pos.add(pos);
                        subTup.superRCs.add(rc);
                        tupList.add(subTup);
                    }
                }
            }
        }
    }
    
    public ArrayList<TreeMap<Integer, T>> getInteractions(){
        return this.interactions;
    }
    public ArrayList<TreeMap<Integer, HigherTupleFinder<T>>> getHigher(){
        return this.higher;
    }
}
