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

package edu.duke.cs.osprey.markstar.visualizer;

import org.ojalgo.matrix.transformation.Rotation;

import java.math.BigDecimal;
import java.util.*;

public class KStarTreeManipulator {

    public static KStarTreeNode consolidateTree(KStarTreeNode root, Collection<String> residues) {
        int[] defaultConf = residues.stream().mapToInt((a)->-1).toArray();
        int[] levelsToKeep =  new int[residues.size()];
        getKeptLevels(0, root, residues, levelsToKeep);
        String[] defaultAssignment = residues.toArray(new String[residues.size()]);
        KStarTreeNode newRoot = new KStarTreeNode(root.level, defaultAssignment, defaultConf,
                root.getLowerBound(), root.getUpperBound(), root.getConfLowerBound(), root.getConfUpperBound(),
                root.epsilon.doubleValue());
        mergeTreeToLevels(newRoot, root, levelsToKeep, 0);
        return newRoot;
    }

    private static void getKeptLevels(int curIndex, KStarTreeNode curNode, Collection<String> residues, int[] keptLevels) {
        String[] assignments = curNode.getAssignments();
        for(int i = 0; i <assignments.length; i++)
        {
            String matched = "";
            for(String residue: residues)
                if(assignments[i].contains(residue) && !assignments[i].contains("*")){
                    keptLevels[curIndex] = curNode.level;
                    matched = residue;
                    curIndex++;
                    break;
                }
            if(!matched.equals(""))
                residues.remove(matched);
        }
        if(curNode.children != null && curNode.children.size() >0)
            getKeptLevels(curIndex, curNode.children.get(0), residues, keptLevels);
    }

    private static void mergeTreeToLevels(KStarTreeNode newRoot, KStarTreeNode oldRoot, int[] levelsToKeep, int curLevelIndex) {
        if(curLevelIndex > levelsToKeep.length)
            return;
        if(oldRoot.children == null || oldRoot.children.size() < 1)
            return;
        int nextLevel = curLevelIndex;
        if(oldRoot.level == levelsToKeep[curLevelIndex])
            nextLevel++;
        for(KStarTreeNode oldChild:oldRoot.children) {
            mergeTreeToLevels(newRoot, oldChild, levelsToKeep, nextLevel);
        }
        if(oldRoot.level == levelsToKeep[curLevelIndex]) {
            Collection<KStarTreeNode> newChildren = mergeSubtreeIntoOneLevel(newRoot, levelsToKeep[curLevelIndex], levelsToKeep);
            newRoot.children.addAll(newChildren);
        }
    }


    public static Collection<KStarTreeNode> mergeChildren(KStarTreeNode subtreeRoot, int targetLevel, int[] indicesToKeep) {

        return null;
    }

    public static Collection<KStarTreeNode> mergeSubtreeIntoOneLevel(KStarTreeNode subtreeRoot, int targetLevel, int[] indicesToKeep) {
        Collection<KStarTreeNode> newChildren = new ArrayList<>();
        Map<Integer, List<KStarTreeNode>> binnedNodes = binNodes(subtreeRoot, targetLevel);
        for(Integer assignmentAtLevel:binnedNodes.keySet())
        {
            int[] newConfAssignments = extractElements(binnedNodes.get(assignmentAtLevel).get(0).getConfAssignments(), indicesToKeep);
            assert(assignmentAtLevel == newConfAssignments[subtreeRoot.level+1]);
            String[] newAssignments = extractElements(binnedNodes.get(assignmentAtLevel).get(0).getAssignments(), indicesToKeep);
            BigDecimal cumulativeLowerBound = binnedNodes.get(assignmentAtLevel).stream().map(KStarTreeNode::getLowerBound).reduce(BigDecimal.ZERO, BigDecimal::add);
            BigDecimal cumulativeUpperBound = binnedNodes.get(assignmentAtLevel).stream().map(KStarTreeNode::getUpperBound).reduce(BigDecimal.ZERO, BigDecimal::add);
            double minConfLower = binnedNodes.get(assignmentAtLevel).stream().map(KStarTreeNode::getConfLowerBound).reduce(Double.POSITIVE_INFINITY, Double::min);
            double maxConfUpper = binnedNodes.get(assignmentAtLevel).stream().map(KStarTreeNode::getConfUpperBound).reduce(Double.NEGATIVE_INFINITY, Double::max);
            KStarTreeNode newNode = new KStarTreeNode(subtreeRoot.level+1, newAssignments, newConfAssignments, cumulativeLowerBound, cumulativeUpperBound,
                    minConfLower, maxConfUpper, subtreeRoot.epsilon.doubleValue());
            newChildren.add(newNode);
        }
        return newChildren;
    }

    private static Map<Integer, List<KStarTreeNode>> binNodes(KStarTreeNode subtreeRoot, int targetLevel) {
        Map<Integer, List<KStarTreeNode>> binnedNodes = new HashMap<>();
        binNodes(subtreeRoot, binnedNodes, targetLevel);
        return binnedNodes;
    }

    public static KStarTreeNode mergeSubtreeLevels(KStarTreeNode subtreeRoot, int targetLevel) {
        List<KStarTreeNode> children = subtreeRoot.children;
        if(children == null || children.size() < 1)
            return null;
        Map<Integer, List<KStarTreeNode>> binnedNodes = new HashMap<>();
        binNodes(subtreeRoot, binnedNodes, targetLevel);
        int skippedLevels = targetLevel - subtreeRoot.level;

        int[] splicedRootConfAssignments = cutConfAssignments(subtreeRoot.getConfAssignments(), subtreeRoot.level, targetLevel-1);
        String[] splicedRootAssignments = cutStringAssignments(subtreeRoot.getAssignments(), subtreeRoot.level, targetLevel-1);
        KStarTreeNode newSubtreeRoot = new KStarTreeNode(subtreeRoot.level, splicedRootAssignments, splicedRootConfAssignments,
                subtreeRoot.getLowerBound(), subtreeRoot.getUpperBound(), subtreeRoot.getConfLowerBound(), subtreeRoot.getConfUpperBound(),
                subtreeRoot.epsilon.doubleValue());
        int[] confAssignments = subtreeRoot.getConfAssignments();
        for(Integer assignmentAtLevel:binnedNodes.keySet())
        {
            int[] splicedConfAssignments = cutConfAssignments(confAssignments, subtreeRoot.level, targetLevel-1);
            String[] splicedAssignments = cutStringAssignments(subtreeRoot.getAssignments(), subtreeRoot.level, targetLevel-1);
            BigDecimal cumulativeLowerBound = binnedNodes.get(assignmentAtLevel).stream().map(KStarTreeNode::getLowerBound).reduce(BigDecimal.ZERO, BigDecimal::add);
            BigDecimal cumulativeUpperBound = binnedNodes.get(assignmentAtLevel).stream().map(KStarTreeNode::getUpperBound).reduce(BigDecimal.ZERO, BigDecimal::add);
            double minConfLower = binnedNodes.get(assignmentAtLevel).stream().map(KStarTreeNode::getConfLowerBound).reduce(Double.POSITIVE_INFINITY, Double::min);
            double maxConfUpper = binnedNodes.get(assignmentAtLevel).stream().map(KStarTreeNode::getConfUpperBound).reduce(Double.NEGATIVE_INFINITY, Double::max);
            KStarTreeNode newNode = new KStarTreeNode(subtreeRoot.level, splicedAssignments, splicedConfAssignments, cumulativeLowerBound, cumulativeUpperBound,
                    minConfLower, maxConfUpper, subtreeRoot.epsilon.doubleValue());
            newSubtreeRoot.addChild(newNode);
        }
        return newSubtreeRoot;
    }


    private static String[] extractElements(String[] source, int[] indices) {
        String[] out = new String[indices.length];
        for (int i = 0; i < indices.length; i++) {
            out[i] = source[indices[i]];
        }
        return out;
    }

    private static int[] extractElements(int[] source, int[] indices) {
        int[] out = new int[indices.length];
        for (int i = 0; i < indices.length; i++) {
            out[i] = source[indices[i]];
        }
        return out;
    }

    private static int[] cutConfAssignments(int[] confAssignments, int firstCutIndex, int lastCutIndex) {
        Integer[] objectArray = new Integer[confAssignments.length];
        return  Arrays.stream((Integer[])replaceRangeWithEmptyArrayElement(
                objectArray, firstCutIndex, lastCutIndex))
                .mapToInt(Integer::intValue)
                .toArray();
    }

    private static String[] cutStringAssignments(String[] assignments, int firstCutIndex, int lastCutIndex) {
        return (String[]) replaceRangeWithEmptyArrayElement(assignments, firstCutIndex, lastCutIndex);
    }

    private static Object[] replaceRangeWithEmptyArrayElement(Object[] source, int startCut, int endCut) {
        Object[] out = new Object[source.length - (endCut - startCut)];
        System.arraycopy(source, 0, out, 0, startCut);
        System.arraycopy(source, endCut+1, out, startCut+2, source.length - endCut);
        return out;
    }

    private static void binNodes(KStarTreeNode subtreeRoot, Map<Integer, List<KStarTreeNode>> binnedNodes, int targetLevel) {
        int level = subtreeRoot.level;
        if(level > targetLevel)
            return;
        if(level == targetLevel)
        {
            int[] confAssignments = subtreeRoot.getConfAssignments();
            int assignmentAtLevel = confAssignments[targetLevel-1];
            if(!binnedNodes.containsKey(assignmentAtLevel))
                binnedNodes.put(assignmentAtLevel, new ArrayList<>());
            binnedNodes.get(assignmentAtLevel).add(subtreeRoot);
            return;
        }

        List<KStarTreeNode> children = subtreeRoot.children;
        if(children == null || children.size() < 1)
            return;
        for(KStarTreeNode child: children)
            binNodes(child, binnedNodes, targetLevel);
    }
}
