package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class TupleTrie<T extends RCTupleContainer> {
    public final static int WILDCARD_RC = -123;
    public static boolean debug = false;
    TupleTrieNode root;
    List<SimpleConfSpace.Position> positions;
    private int numEntries;

    public TupleTrie(List<SimpleConfSpace.Position> positions)
    {
        this.positions = positions;
        this.root = new TupleTrieNode(positions, -1);
    }

    public void insert(T entry) {
        if(debug)
            checkRCTuple(entry.tup);
        root.insert(entry, 0);
        numEntries++;
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

    public List<T> getEntries(RCTuple query) {
        List<T> entries = new ArrayList<>();
        root.populateEntries(query.sorted(), entries);
        return entries;
    }

    public boolean contains(RCTuple query) {
        return root.contains(query.sorted(), 0);

    }

    public int size() {
        return numEntries;
    }

    public List<T> getAllEntries(){
        List<T> output = new ArrayList<>();
        root.getAllEntries(output);
        return output;
    }

    public void clear() {
        this.root = new TupleTrieNode(positions, -1);
    }

    public void writeEntriesToFile(String filename){
        File file = new File(filename);
        writeEntriesToFile(file);
    }
    public void writeEntriesToFile(File f) {
        List<String> lineList = getAllEntries().stream()
                .distinct().map(T::toString).collect(Collectors.toList());
        FileTools.writeFile(String.join("\n", lineList), f);
    }

    private class TupleTrieNode {
        // Wildcard rc
        int rc = WILDCARD_RC;
        int positionIndex = -1;
        int position = -1;
        List<SimpleConfSpace.Position> positions;
        List<T> entries = new ArrayList<>();
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

        public void insert(T entry, int tupIndex)
        {
            for(TupleTrieNode child: children.values()) {
                debugPrint(this+"->"+child);
            }
            for(T e: entries)
            {
                debugPrint(e.toString());
            }
            RCTuple tup = entry.tup;
            if(tupIndex >= tup.size()) {
                debugPrint("Reached end of tuple, inserting correction at "+this+".");
                entries.add(entry);
                for(T e: entries)
                {
                    debugPrint(e.toString());
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
                children.get(WILDCARD_RC).insert(entry, tupIndex);
            }
            else
            {
                if(!children.containsKey(childRC)) {
                    TupleTrieNode newChild = new TupleTrieNode(positions, positionIndex+1);
                    newChild.rc = childRC;
                    children.put(childRC, newChild);
                    debugPrint("Added child "+newChild+" to "+this);
                }
                children.get(childRC).insert(entry, tupIndex+1);
            }

        }


        public void populateEntries (RCTuple query, List<T> output) {
            debugPrint("Matching entries for "+query.stringListing());
            populateEntries(query, output, 0);
        }

        private void populateEntries(RCTuple query, List<T> output, int tupleIndex) {
            debugPrint("Currently at "+this);
            if(entries.size() > 0)
            {
                output.addAll(entries);
                debugPrint("Adding entries from "+this);
                if(debug)
                    for(T entry:entries) {
                        System.out.println(entry);
                    }
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
            if(position + 1 ==  currentPos && children.containsKey(currentRC))
                children.get(currentRC).populateEntries(query, output, nextIndex);
            // Also branch on wildcard.
            if(!children.containsKey(WILDCARD_RC))
                children.put(WILDCARD_RC, new TupleTrieNode(positions, positionIndex+1));
            children.get(WILDCARD_RC).populateEntries(query, output, nextIndex);
        }

        public void getAllEntries(List<T> output){
            if(entries.size() > 0)
            {
                output.addAll(entries);
                debugPrint("Adding entries from "+this);
            }

            for (TupleTrieNode child : this.children.values()){
                child.getAllEntries(output);
            }
        }

    }
}
