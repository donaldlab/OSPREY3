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

    public static List<Map<String,List<BigDecimal>>> calcResidueOccupancyList(KStarTreeNode rootNode){
        /**
         * Calculates the occupancy of each rotamer for all residues
         */
        BigDecimal overallUpperBound = rootNode.getUpperBound();
        BigDecimal overallLowerBound = rootNode.getLowerBound();

        List<Map<String,List<BigDecimal>>> marginTree = marginalizeTree(rootNode);

        List<Map<String, List<BigDecimal>>> occTree = new ArrayList<>();

        for( Map<String, List<BigDecimal>> residue: marginTree ) {
            occTree.add(residue.entrySet()
                    .stream().collect(Collectors.toMap(Map.Entry::getKey,
                            e -> Arrays.asList(
                                    e.getValue().get(0).divide(overallUpperBound,10, RoundingMode.HALF_UP),
                                    e.getValue().get(1).divide(overallLowerBound,10, RoundingMode.HALF_UP)
                            ))));
        }

        return occTree;
    }
    public static void printOccupancyList(List<Map<String,List<BigDecimal>>> occTree){
        for( Map<String,List<BigDecimal>> residue : occTree){
            System.out.println(String.format("Residue %d has lower bounded occupancies:\n\t%s",
                    occTree.indexOf(residue),
                    residue.toString()));
        }
        for( Map<String,List<BigDecimal>> residue : occTree){
            System.out.println(String.format("Residue %d has upper bounded occupancies:\n\t%s",
                    occTree.indexOf(residue),
                    residue.toString()));
        }
    }
    public static List<Map<String,List<BigDecimal>>> marginalizeTree(KStarTreeNode rootNode){
        /**
         * Marginalizes a tree by residue
         *
         * This method takes a KStar tree (representing a Boltzmann distribution) and marginalizes
         * the distribution over each residue.
         *
         * TODO: Fix issue where sum of marginal pfuncs exceeds overall pFunc bounds
         * I think this is something to do with rounding, but I believe it is possible that it may be more complicated
         */

        List<Map<String,List<BigDecimal>>> marginalizedTree = new ArrayList<>();

        List<Map<String,List<KStarTreeNode>>> flatTree = KStarTreeManipulator.binAllNodes(rootNode);

        for( Map<String, List<KStarTreeNode>> level : flatTree){
            Map<String, List<BigDecimal>> newResidue = new HashMap<>();
            // for each rotamer at a residue, accumulate the upper and lower bounds
            for( String rotamer : level.keySet() ){
                //generate a pair (Lower, Upper) of bounds for each marginalized rotamer
                newResidue.put(
                        rotamer,Arrays.asList(
                            level.get(rotamer).stream()
                                .map(KStarTreeNode::getLowerBound)
                                .reduce(BigDecimal.ZERO, BigDecimal::add),
                            level.get(rotamer).stream()
                                    .map(KStarTreeNode::getUpperBound)
                                    .reduce(BigDecimal.ZERO, BigDecimal::add)));

            }
            // each map contains upper and lower bounds marginalized for a residue.
            marginalizedTree.add(newResidue);
        }
        // test to make sure we aren't missing significant amounts of pfunc
        testMarginalizedTree(marginalizedTree,rootNode.getUpperBound(),rootNode.getLowerBound());
        return marginalizedTree;
    }
    public static void testMarginalizedTree(List<Map<String,List<BigDecimal>>> mTree, BigDecimal upperBound, BigDecimal lowerBound, Boolean quiet){
        /**
         * Tests to ensure that each marginal distribution contained in a marginalized Tree captures the full Boltzmann distribution
         *
         * Something has likely gone wrong if this is not the case (although it is possible that we have not explored nodes for some residues
         * Also I THINK that we have rounding errors because of the necessity of printing things out to file
         */

        for( Map<String, List<BigDecimal>> residue : mTree) {
            List<BigDecimal> totalStatWeight= residue.values().stream()
                    .reduce(Arrays.asList(BigDecimal.ZERO,BigDecimal.ZERO),
                            (a, b) -> Arrays.asList(
                                    a.get(0).add(b.get(0)),
                                    a.get(1).add(b.get(1))
                            ));
            if (!quiet) {
                System.out.println(String.format("Residue %d has marginal bounds of [%.8e, %.8e](L,U). Full pfunc -> [%.8e, %.8e](L,U)",
                        mTree.indexOf(residue),
                        totalStatWeight.get(0),
                        totalStatWeight.get(1),
                        lowerBound,
                        upperBound
                ));
            }
            List<Double> percentDiff = Arrays.asList(
                    totalStatWeight.get(0).subtract(lowerBound).divide(lowerBound, 10, BigDecimal.ROUND_HALF_UP ).doubleValue(),
                    totalStatWeight.get(1).subtract(upperBound).divide(upperBound, 10, BigDecimal.ROUND_HALF_UP ).doubleValue()
                    );
            if (Math.abs(percentDiff.get(0)) >= 0.000001 || Math.abs(percentDiff.get(1)) >= 0.000001){
                System.err.println(String.format("WARNING: Residue %d marginals are off of the full pfunc by [%.8f, %.8f](L,U) %%",
                        mTree.indexOf(residue),
                        percentDiff.get(0)*100,
                        percentDiff.get(1)*100
                        ));
            }

        }
    }

    public static void testMarginalizedTree(List<Map<String,List<BigDecimal>>> mTree, BigDecimal upperBound, BigDecimal lowerBound){
        testMarginalizedTree(mTree, upperBound, lowerBound, true);
    }
}
