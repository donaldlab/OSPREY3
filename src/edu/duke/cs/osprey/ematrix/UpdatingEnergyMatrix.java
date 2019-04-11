/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.confspace.TupE;

import java.util.*;


public class UpdatingEnergyMatrix extends ProxyEnergyMatrix {
    // Store the seen confs in a trie with wildcards.
    private static final boolean debug = false;
    private TupleTrie corrections;
    private int numPos;
    
    //debug variable
    public final ConfEnergyCalculator sourceECalc;

    public UpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target, ConfEnergyCalculator confECalc) {
        super(confSpace, target);
        corrections = new TupleTrie(confSpace.positions);
        this.numPos = confSpace.getNumPos();
        this.sourceECalc = confECalc;

    }

    public UpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
        super(confSpace, target);
        this.numPos = confSpace.getNumPos();
        this.sourceECalc = null;
        corrections = new TupleTrie(confSpace.positions);
    }

    /*Hack 1: Don't share residues*/
    @Override
    public boolean hasHigherOrderTerms() {
        return corrections.size() > 0;
    }

    @Override
    public boolean hasHigherOrderTermFor(RCTuple query) {
        return corrections.contains(query);
    }

    public String formatCorrections(List<TupE> corrections) {
        String out = "";
        for(TupE correction:corrections)
            out+=correction.tup.stringListing()+":"+correction.E+"\n";
            return out;
    }


    public double getInternalEnergy(RCTuple tup){
        //internal energy of a tuple of residues when they're in the specified RCs

        // OPTIMIZATION: don't even check higher terms if the energy matrix doesn't have any
        // this does wonders to CPU cache performance!
        boolean useHigherOrderTerms = hasHigherOrderTerms();

        ArrayList<Integer> tuppos = tup.pos;
        ArrayList<Integer> tupRCs = tup.RCs;

        // OPTIMIZATION: split oneBody and pairwise energies into separate loops
        // to improve CPU cache performance

        int numPosInTuple = tup.pos.size();
        double energy = 0;

        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            int posNum = tuppos.get(indexInTuple);
            int RCNum = tupRCs.get(indexInTuple);

            energy += getOneBody(posNum,RCNum);
        }

        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            int posNum = tuppos.get(indexInTuple);
            int RCNum = tupRCs.get(indexInTuple);

            for(int index2=0; index2<indexInTuple; index2++){
                int pos2 = tuppos.get(index2);
                int rc2 = tupRCs.get(index2);

                energy += getPairwise(posNum,RCNum,pos2,rc2);

            }
        }
        if (useHigherOrderTerms) {
            energy+=internalEHigherOrder(tup);
        }

        return energy;
    }


    /*
    public double getPairwiseE(int res1, int AA1, int rot1, int res2, int AA2, int rot2){
    //lookup by residue number, amino-acid index (for the residue's rotamer library),
    //and RC index for that residue and AA type (rotamer index for AA type if rigid backbone,
    //otherwise defined in the ConfSearchSpace)
    //RETURN ERROR IF RES1 AND RES2 ARE NOT SINGLE RESIDUES
    return getPairwise(res1, index1, res2, index2);
    }
     */
    //intra+shell similar...

    double internalEHigherOrder(RCTuple tup){
        //Computes the portion of the internal energy for tuple tup
        //that consists of interactions in htf (corresponds to some sub-tuple of tup)
        //with RCs whose indices in tup are < curIndex
        double E = 0;
        //E += super.internalEHigherOrder(tup, curIndex, htf);
        List<TupE> confCorrections = corrections.getCorrections(tup);
        if(confCorrections.size() > 0) {
            double corr = processCorrections(confCorrections);
            E += corr;
        }

        return E;
    }

    private double processCorrections(List<TupE> confCorrections) {
        Collections.sort(confCorrections, (a,b)->-Double.compare(a.E,b.E));
        double sum = 0;
        // Attempt 1: be greedy and start from the largest correction you
        // can get instead of trying to solve the NP-Complete problem.
        Set<Integer> usedPositions = new HashSet<>();
        List<TupE> usedCorrections = new ArrayList<>();
        int numApplied = 0;
        for(TupE correction: confCorrections) {
            if (usedPositions.size() >= numPos) {
                break;
            }
            Collection<Integer> positions = correction.tup.pos;
            boolean noIntersections = true;
            for(int position : positions) {
                if(usedPositions.contains(position)) {
                    noIntersections = false;
                    break;
                }
            }
            if(noIntersections) {
                usedPositions.addAll(correction.tup.pos);
                usedCorrections.add(correction);
                //System.out.println("Applying correction "+correction.tup.stringListing()+":"+correction.E);
                numApplied++;
                sum += correction.E;
            }
        }
        if(debug)
        {
            Set<Integer> positionCheck = new HashSet<>();
            for(TupE correction: usedCorrections) {
                for(int pos: correction.tup.pos) {
                    if (positionCheck.contains(pos))
                        System.err.println("REUSING POSITION "+pos);
                }
            }
        }
        return sum;
    }

    @Override
    public void setHigherOrder(RCTuple tup, Double val) {
        if(tup.size() < 3)
        {
            System.err.println("Should not be trying to submit correction of lower-order term.");
            return;
        }
        RCTuple orderedTup = tup.sorted();
        corrections.insert(new TupE(orderedTup, val));
    }

    public static class TupleTrie {
        public final static int WILDCARD_RC = -123;
        TupleTrieNode root;
        List<SimpleConfSpace.Position> positions;
        private int numCorrections;
        public TupleTrie(List<SimpleConfSpace.Position> positions)
        {
            this.positions = positions;
            root = createTrie(positions);
        }

        private TupleTrieNode createTrie(List<SimpleConfSpace.Position> positions) {
            root = new TupleTrieNode(positions, -1);
            return root;
        }

        public void insert(TupE correction) {
            if(debug)
                checkRCTuple(correction.tup);
            root.insert(correction, 0);
            numCorrections++;
        }

        private void checkRCTuple(RCTuple tup) {
            int lastIndex = 0;
            for(int i = 0; i < tup.size(); i++)
            {
                int index = positions.indexOf(tup.pos.get(i));
                if(index > -1 && lastIndex > index)
                    System.err.println("Tuple and confspace are not ordered the same way.");
                lastIndex = index;
            }
        }

        public List<TupE> getCorrections(RCTuple query) {
            List<TupE> corrections = new ArrayList<>();
            root.populateCorrections(query.sorted(), corrections);
            return corrections;
        }

        public boolean contains(RCTuple query) {
            return root.contains(query.sorted(), 0);

        }

        public int size() {
            return numCorrections;
        }


        private class TupleTrieNode {
            // Wildcard rc
            int rc = WILDCARD_RC;
            int positionIndex = -1;
            int position = -1;
            List<SimpleConfSpace.Position> positions;
            List<TupE> corrections = new ArrayList<>();
            Map<Integer, TupleTrieNode> children = new HashMap<>();

            private TupleTrieNode(List<SimpleConfSpace.Position> positions, int positionIndex) {
                this.positions = positions;
                this.positionIndex = positionIndex;
                if(positionIndex >= 0)
                    this.position = positions.get(positionIndex).index;
                if(positionIndex+1 < positions.size())
                    children.put(WILDCARD_RC, new TupleTrieNode(positions, positionIndex+1));
            }

            public boolean contains(RCTuple query, int tupleIndex) {
                /*
                if(query.size() != positions.size())
                    System.err.println("Querying corrections for a partial conf. This is likely unintentional.");
                    */
                debugPrint("Currently at "+this);
                if(tupleIndex >= query.size())
                    return true;
                int currentRC = query.RCs.get(tupleIndex);
                int currentPos = query.pos.get(tupleIndex);
                int indexedPos = -1;
                int indexedRC = WILDCARD_RC;
                if(tupleIndex > 0) {
                    indexedRC = query.RCs.get(tupleIndex-1);
                    indexedPos = query.pos.get(tupleIndex-1);
                }
                if(tupleIndex + 1 == positions.size())
                    return true;
                int nextIndex = tupleIndex + 1;
                if(position + 1 < currentPos) {
                    if(children.get(WILDCARD_RC) == null)
                        return false;
                    return children.get(WILDCARD_RC).contains(query, tupleIndex);
                }
                if(!children.containsKey(currentRC))
                    return false;
                if(children.get(currentRC) == null)
                    return false;
                return children.get(currentRC).contains(query, nextIndex);
            }

            public String toString()
            {
                String rcString = "*";
                if(rc > -1)
                    rcString = ""+rc;
                return ""+position+":"+rcString;
            }

            private void debugPrint(String s)
            {if(debug)System.out.println(s);}

            public void insert(TupE correction, int tupIndex)
            {
                for(TupleTrieNode child: children.values()) {
                    debugPrint(this+"->"+child);
                }
                for(TupE corr: corrections)
                {
                   debugPrint(corr.tup.stringListing()+":"+corr.E);
                }
                RCTuple tup = correction.tup;
                if(tupIndex >= tup.size()) {
                    debugPrint("Reached end of tuple, inserting correction at "+this+".");
                    corrections.add(correction);
                    for(TupE corr: corrections)
                    {
                        debugPrint(corr.tup.stringListing()+":"+corr.E);
                    }
                    return;
                }
                int currentIndex = -1;
                int nodeIndex = position;
                int currentRC = WILDCARD_RC;
                if(tupIndex > 0) {
                    currentIndex = tup.pos.get(tupIndex - 1);
                    currentRC = tup.pos.get(tupIndex - 1);
                }
                int childIndex = tup.pos.get(tupIndex);
                int childRC = tup.RCs.get(tupIndex);
                if(nodeIndex+1 != childIndex) {
                   debugPrint((nodeIndex+1)+"!="+childIndex+", continuing...");
                    children.get(WILDCARD_RC).insert(correction, tupIndex);
                }
                else
                {
                    if(!children.containsKey(childRC)) {
                        TupleTrieNode newChild = new TupleTrieNode(positions, positionIndex+1);
                        newChild.rc = childRC;
                        children.put(childRC, newChild);
                        debugPrint("Added child "+newChild+" to "+this);
                    }
                    children.get(childRC).insert(correction, tupIndex+1);
                }

            }


            public void populateCorrections (RCTuple query, List<TupE> output) {
                debugPrint("Matching corrections for "+query.stringListing());
                populateCorrections(query, output, 0);
            }

            private void populateCorrections(RCTuple query, List<TupE> output, int tupleIndex) {
                debugPrint("Currently at "+this);
                if(corrections.size() > 0)
                {
                    output.addAll(corrections);
                    debugPrint("Adding corrections from "+this);
                }
                if(tupleIndex >= query.size())
                    return;
                int currentRC = query.RCs.get(tupleIndex);
                int currentPos = query.pos.get(tupleIndex);
                int indexedPos = -1;
                int indexedRC = WILDCARD_RC;
                if(tupleIndex > 0) {
                    indexedRC = query.RCs.get(tupleIndex-1);
                    indexedPos = query.pos.get(tupleIndex-1);
                }
                if(indexedPos > position || (indexedPos == position && indexedRC!= rc && rc != WILDCARD_RC))
                    System.err.println("Error in trie traversal.");
                if(tupleIndex + 1 > positions.size())
                    return;
                int nextIndex = tupleIndex + 1;
                if(position + 1 < currentPos)
                    nextIndex = tupleIndex;
                if(children.containsKey(currentRC))
                    children.get(currentRC).populateCorrections(query, output, nextIndex);
                // Also branch on wildcard.
                if(!children.containsKey(WILDCARD_RC))
                    children.put(WILDCARD_RC, new TupleTrieNode(positions, positionIndex+1));
                children.get(WILDCARD_RC).populateCorrections(query, output, nextIndex);
            }
        }

    }

}
