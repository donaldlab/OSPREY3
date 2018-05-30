package edu.duke.cs.osprey.astar.ewakstar;

import java.util.*;


public class EWAKStarLimitedSequenceTrie {

    private static final TrieNode[] EMPTYNODES = new TrieNode[0];

    private static final class TrieNode {

        private final String aminoAcid;
        private boolean isSeq = false;
        private Map<String, TrieNode> children = null;

        public TrieNode(String aa) {
            aminoAcid = aa;
        }

        public TrieNode getOrCreateChild(String aa) {
            if (children == null) {
                children = new HashMap<>();
            }
            TrieNode kid = children.get(aa);
            if (kid == null) {
                kid = new TrieNode(aa);
                children.put(aa, kid);
            }

            return kid;
        }

        public TrieNode get(String aa) {
            return children != null ? children.get(aa) : null;
        }

        public void setSeq() {
            isSeq = true;
        }

        public boolean isSeq() {
            return isSeq;
        }

        public String getAA() {
            return aminoAcid;
        }

        public TrieNode[] getChildNodes() {
            if (children == null) {
                return EMPTYNODES;
            }
            TrieNode[] result = children.values().toArray(new TrieNode[children.size()]);
            return result;
        }

    }

    private final TrieNode root;
    private int size = 0;

    public EWAKStarLimitedSequenceTrie(){
        // root has null character.
        root = new TrieNode("null");
    }

    public void addSeq(String sequence){
        TrieNode node = root;
        for (String aa : sequence.split(" ")) {
            node = node.getOrCreateChild(aa);
        }
        if (!node.isSeq()) { // fix - only add new words....
            node.setSeq();
            size++;
        }
    }

    public void printSize(){
        System.out.println(size);
    }
    public boolean containsSeq(String sequence){
        TrieNode node = root;
        for (String aa : sequence.split(" ")) {
            node = node.get(aa);
            if (node == null) {
                break;
            }
        }
        return node != null;
    }

    public Set<String> getSeq(String sequence){
        TrieNode node = root;
        for (String aa : sequence.split(" ")) {
            node = node.get(aa);
        }
        Set<String> resTypes = new HashSet<>();
        for(TrieNode n:node.getChildNodes()){
            resTypes.add(n.getAA());
        }
        return resTypes;
    }
}