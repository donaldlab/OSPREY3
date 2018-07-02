package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupE;

import java.util.*;


public class UpdatingEnergyMatrix extends ProxyEnergyMatrix {
    // Store the seen confs in a trie with wildcards.
    private static final boolean debug = false;
    private TupleTrie corrections;
    private int numPos;

    public UpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
        super(confSpace, target);
        this.numPos = confSpace.getNumPos();
        corrections = new TupleTrie(confSpace.positions);
    }

    /*Hack 1: Don't share residues*/

    @Override
    double internalEHigherOrder(RCTuple tup, int curIndex, HigherTupleFinder<Double> htf){
        //Computes the portion of the internal energy for tuple tup
        //that consists of interactions in htf (corresponds to some sub-tuple of tup)
        //with RCs whose indices in tup are < curIndex
        double E = 0;
        E += super.internalEHigherOrder(tup, curIndex, htf);

        List<TupE> confCorrections = corrections.getCorrections(tup);
        E += processCorrections(confCorrections);

        return E;
    }

    private double processCorrections(List<TupE> confCorrections) {
        Collections.sort(confCorrections);
        double sum = 0;
        // Attempt 1: be greedy and start from the largest correction you
        // can get instead of trying to solve the NP-Complete problem.
        Set<Integer> usedPositions = new HashSet<>();
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
                sum += correction.E;
            }
        }
        return sum;
    }

    @Override
    public void setHigherOrder(RCTuple tup, Double val) {
        super.setHigherOrder(tup, val);
        RCTuple orderedTup = tup.orderByPos();
        corrections.insert(new TupE(orderedTup, val));
    }

    public void addMinimizedConf(RCTuple tup) {

    }

    public static class TupleTrie {
        public final static int WILDCARD_RC = -123;
        TupleTrieNode root;
        List<SimpleConfSpace.Position> positions;
        public TupleTrie(List<SimpleConfSpace.Position> positions)
        {
            this.positions = positions;
            root = createTrie(positions);
        }

        private TupleTrieNode createTrie(List<SimpleConfSpace.Position> positions) {
            root = new TupleTrieNode(positions, -1);
            TupleTrieNode next = root;

            return root;
        }

        public void insert(TupE correction) {
            if(debug)
                checkRCTuple(correction.tup);
            root.insert(correction, 0);
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
            root.populateCorrections(query.orderByPos(), corrections);
            return corrections;
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
                        System.out.println(corr.tup.stringListing()+":"+corr.E);
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
                /*
                if(query.size() != positions.size())
                    System.err.println("Querying corrections for a partial conf. This is likely unintentional.");
                    */
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
                if(tupleIndex + 1 == positions.size())
                    return;
                int nextIndex = tupleIndex + 1;
                if(position + 1 < currentPos)
                    nextIndex = tupleIndex;
                if(children.containsKey(currentRC))
                    children.get(currentRC).populateCorrections(query, output, nextIndex);
                // Also branch on wildcard.
                children.get(WILDCARD_RC).populateCorrections(query, output, nextIndex);
            }
        }

    }

}
