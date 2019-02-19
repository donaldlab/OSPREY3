package edu.duke.cs.osprey.markstar.visualizer;

import java.math.BigDecimal;
import java.util.*;

public class KStarTreeAnalyzer {

    public static void calcResidueEntropy(KStarTreeNode rootNode){
        /**
         * Calculates residue entropy by collecting nodes for each residue, treating this node as root
         */
        //call something that calculates bounds on residue occupancy
        //do calculations to turn this into entropy bounds
    }

    public static void calcResidueOccupancy(KStarTreeNode rootNode){
        /**
         * Calculates the occupancy of each rotamer for all residues
         */
        BigDecimal overallUpperBound = rootNode.getUpperBound();
        BigDecimal overallLowerBound = rootNode.getLowerBound();
        List<Map<Integer, List<KStarTreeNode>>> flatTree = KStarTreeManipulator.binAllNodes(rootNode);
        List<Map<Integer, List<KStarTreeNode>>> reducedTree = new ArrayList<>();

        for( Map<Integer, List<KStarTreeNode>> level : flatTree){
            for( Map.Entry<Integer,List<KStarTreeNode>> assList : level.entrySet() ){
                assList.getValue().stream();//TODO: map and reduce
            }
        }

    }
}
