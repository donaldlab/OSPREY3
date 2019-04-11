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

package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.confspace.SeqSpace;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;


public class EwakstarLimitedSequenceTrie {

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
    public final SeqSpace seqSpace;
    private int size = 0;

    public EwakstarLimitedSequenceTrie(SeqSpace ss){
        // root has null character.
        root = new TrieNode("null");
        seqSpace = ss;
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

    public Set<String> getFirstPos(){
        Set<String> resTypes = new HashSet<>();
        for(TrieNode n: root.getChildNodes()){
            resTypes.add(n.getAA().split("=")[1].toUpperCase());
        }

        return resTypes;
    }

    public Set<String> getSeq(String sequence){
        TrieNode node = root;
        for (String aa : sequence.split(" ")) {
            node = node.get(aa);
        }
        Set<String> resTypes = new HashSet<>();
        for(TrieNode n:node.getChildNodes()){
            resTypes.add(n.getAA().split("=")[1].toUpperCase());
        }
        return resTypes;
    }
}
