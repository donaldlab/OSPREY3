package edu.duke.cs.osprey.markstar.visualizer;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;
import java.util.stream.Collectors;

public class KStarTreeAnalyzer {

    public static void calcResidueEntropy(KStarTreeNode rootNode){
        /**
         * Calculates residue entropy by collecting nodes for each residue, treating this node as root
         */
        //call something that calculates bounds on residue occupancy
        //do calculations to turn this into entropy bounds
    }

    public static List<List<Map<Integer,BigDecimal>>> calcResidueOccupancyList(KStarTreeNode rootNode){
        /**
         * Calculates the occupancy of each rotamer for all residues
         */
        BigDecimal overallUpperBound = rootNode.getUpperBound();
        BigDecimal overallLowerBound = rootNode.getLowerBound();

        List<List<Map<Integer,BigDecimal>>> marginTree = marginalizeTree(rootNode);

        List<List<Map<Integer, BigDecimal>>> occTree = new ArrayList<>();

        List<Map<Integer, BigDecimal>> lowerList = new ArrayList<>();
        List<Map<Integer, BigDecimal>> upperList = new ArrayList<>();

        for( Map<Integer, BigDecimal> residueLower : marginTree.get(0) ) {
            lowerList.add(residueLower.entrySet()
                    .stream().collect(Collectors.toMap(Map.Entry::getKey, e -> e.getValue().divide(overallLowerBound,10, RoundingMode.HALF_UP))));
        }
        for( Map<Integer, BigDecimal> residueUpper : marginTree.get(1) ) {
            upperList.add(residueUpper.entrySet()
                    .stream().collect(Collectors.toMap(Map.Entry::getKey, e -> e.getValue().divide(overallUpperBound,10, RoundingMode.HALF_UP))));
        }

        occTree.add(lowerList);
        occTree.add(upperList);
        return occTree;
    }
    public static void printOccupancyList(List<List<Map<Integer,BigDecimal>>> occTree){
        for( Map<Integer,BigDecimal> residue : occTree.get(0)){
            System.out.println(String.format("Residue %d has lower bounded occupancies:\n\t%s",
                    occTree.get(0).indexOf(residue),
                    residue.toString()));
        }
        for( Map<Integer,BigDecimal> residue : occTree.get(1)){
            System.out.println(String.format("Residue %d has upper bounded occupancies:\n\t%s",
                    occTree.get(1).indexOf(residue),
                    residue.toString()));
        }
    }
    public static List<List<Map<Integer,BigDecimal>>> marginalizeTree(KStarTreeNode rootNode){
        /**
         * Marginalizes a tree by residue
         *
         * This method takes a KStar tree (representing a Boltzmann distribution) and marginalizes
         * the distribution over each residue.
         */

        List<List<Map<Integer, BigDecimal>>> boundsList = new ArrayList<>();
        List<Map<Integer, BigDecimal>> reducedTreeLower = new ArrayList<>();
        List<Map<Integer, BigDecimal>> reducedTreeUpper = new ArrayList<>();

        List<Map<Integer, List<KStarTreeNode>>> flatTree = KStarTreeManipulator.binAllNodes(rootNode);

        for( Map<Integer, List<KStarTreeNode>> level : flatTree){
            Map<Integer, BigDecimal> newMapLower = new HashMap<>();
            Map<Integer, BigDecimal> newMapUpper = new HashMap<>();
            // for each rotamer at a residue, accumulate the upper and lower bounds
            for( Integer rotamer : level.keySet() ){
                newMapLower.put(rotamer,
                        level.get(rotamer).stream()
                                .map(KStarTreeNode::getLowerBound)
                                .reduce(BigDecimal.ZERO, BigDecimal::add));
                newMapUpper.put(rotamer,
                        level.get(rotamer).stream()
                                .map(KStarTreeNode::getUpperBound)
                                .reduce(BigDecimal.ZERO, BigDecimal::add));
            }
            // each map contains upper and lower bounds marginalized for a residue.
            // add them to the residue list
            reducedTreeLower.add(newMapLower);
            reducedTreeUpper.add(newMapUpper);
        }
        boundsList.add(reducedTreeLower);
        boundsList.add(reducedTreeUpper);
        return boundsList;
    }
    public static void testMarginalizedTree(List<List<Map<Integer,BigDecimal>>> mTree, BigDecimal upperBound, BigDecimal lowerBound, Boolean quiet){
        /**
         * Tests to ensure that each marginal distribution contained in a marginalized Tree captures the full Boltzmann distribution
         *
         * Something has likely gone wrong if this is not the case (although it is possible that we have not explored nodes for some residues
         * Also I THINK that we have rounding errors because of the necessity of printing things out to file
         */

        for( Map<Integer, BigDecimal> residueLower : mTree.get(0) ) {
            BigDecimal totalStatWeightLower = residueLower.values().stream()
                    .reduce(BigDecimal.ZERO, BigDecimal::add);
            if (!quiet) {
                System.out.println(String.format("Residue %d has total lower bound weight of %.8e, while the total lower bound is %.8e.",
                        mTree.get(0).indexOf(residueLower),
                        totalStatWeightLower,
                        lowerBound));
            }
        }
        for( Map<Integer, BigDecimal> residueUpper : mTree.get(1) ){
            BigDecimal totalStatWeightUpper = residueUpper.values().stream()
                    .reduce(BigDecimal.ZERO, BigDecimal::add);
            if(!quiet){
                System.out.println(String.format("Residue %d has total upper bound weight of %.8e, while the total upper bound is %.8e.",
                        mTree.get(1).indexOf(residueUpper),
                        totalStatWeightUpper,
                        upperBound));
            }
        }

    }

    public static void testMarginalizedTree(List<List<Map<Integer,BigDecimal>>> mTree, BigDecimal upperBound, BigDecimal lowerBound){
        testMarginalizedTree(mTree, upperBound, lowerBound, true);
    }
}
