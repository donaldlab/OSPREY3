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

    public static Map<String,Map<String,List<Double>>> calcResidueOccupancyList(KStarTreeNode rootNode){
        /**
         * Calculates the occupancy of each rotamer for all residues
         */
        BigDecimal overallUpperBound = rootNode.getUpperBound();
        BigDecimal overallLowerBound = rootNode.getLowerBound();

        Map<String,Map<String,List<BigDecimal>>> marginTree = marginalizeTree(rootNode);

        Map<String,Map<String, List<Double>>> occTree = new HashMap<>();

        for( String residue : marginTree.keySet() ) {
            Map<String,List<Double>> residueOccupancy = new HashMap<>();
            Map<String,List<BigDecimal>> residueMarginalSum = getCumulativeMarginalRes(marginTree.get(residue), residue, overallLowerBound, overallUpperBound, false);
            // Note that we are forcing occupancies to always sum to 1 despite rounding errors
            // If we have an undershoot, then I believe this should effectively spread the weight equally over rotamers, which is what we want
            // If we have an overshoot, I think that the weight is shifted slightly more toward high occupancy rotamers, as those have been rounded more aggressively

            for(String rotamer : marginTree.get(residue).keySet()){
                residueOccupancy.put(rotamer,
                        Arrays.asList(
                                // Lower bound on occupancy is the lower bound on stat weight over the upper bound on cumulative stat weight (with this rotamer lower bound)
                                marginTree.get(residue).get(rotamer).get(0).divide(residueMarginalSum.get(rotamer).get(1),10,RoundingMode.HALF_UP).doubleValue(),
                                // Upper bound on occupancy is the upper bound on stat weight over the lower bound on cumulative stat weight (with this rotamer upper bound)
                                marginTree.get(residue).get(rotamer).get(1).divide(residueMarginalSum.get(rotamer).get(0),10,RoundingMode.HALF_UP).doubleValue()));
            }
            occTree.put(residue, residueOccupancy);
        }
        return occTree;
    }
    public static void printOccupancyList(Map<String,Map<String,List<Double>>> occTree){
        for( String residue : occTree.keySet()){
            System.out.println(String.format("Residue %s has bounded occupancies:\n\t%s",
                    residue,
                    occTree.get(residue).toString()));
        }
    }
    public static Map<String,Map<String,List<BigDecimal>>> marginalizeTree(KStarTreeNode rootNode){
        /**
         * Marginalizes a tree by residue
         *
         * This method takes a KStar tree (representing a Boltzmann distribution) and marginalizes
         * the distribution over each residue.
         *
         * TODO: Fix issue where sum of marginal pfuncs exceeds overall pFunc bounds
         * I think this is something to do with rounding, but I believe it is possible that it may be more complicated
         */

        Map<String,Map<String,List<BigDecimal>>> marginalizedTree = new HashMap<>();

        List<Map<String,List<KStarTreeNode>>> flatTree = KStarTreeManipulator.binAllNodes(rootNode);

        for( Map<String, List<KStarTreeNode>> level : flatTree){
            Map<String, List<BigDecimal>> newResidue = new HashMap<>();

            // find the key for this residue - that is - check all nodes and make sure they correspond to the same residue
            String resname = "";
            List<String> presentResidues = level.values().stream()
                    .flatMap(c -> c.stream()) //flatten stream
                    .map(KStarTreeNode::getMargResidue) //map to residue name
                    .distinct() //remdup
                    .collect(Collectors.toList());
            if(presentResidues.size() != 1) {
                System.err.println("We have too many residues!!");
            }else{
                resname = presentResidues.get(0);
            }


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
            marginalizedTree.put(resname, newResidue);
        }
        return marginalizedTree;
    }
    public static Map<String, List<BigDecimal>> getCumulativeMarginalRes(Map<String,List<BigDecimal>> residue, String resname, BigDecimal lowerBound, BigDecimal upperBound, Boolean quiet) {
        /**
         * Tests to ensure that the marginal distribution over a residue captures the full Boltzmann distribution
         *
         * Something has likely gone wrong if this is not the case (although it is possible that we have not explored nodes for some residues
         * Also I THINK that we have rounding errors because of the necessity of printing things out to file
         *
         * Returns the sum of the marginal partition function for a residue
         */

        Map<String, List<BigDecimal>> rotamerCumulativeWeights = new HashMap<>();
        for (String rotamer : residue.keySet()) {
            List<BigDecimal> boundedStatWeight = residue.entrySet().stream()
                    .map(
                            a -> {
                                if (a.getKey().equals(rotamer)) {
                                    // If a is the rotamer for which we are calculating max relative stat weight,
                                    // then the lower bound is equal to the highest of a + the lowest of all others.
                                    // and the upper bound is equal to the lowest of a + the highest of all others,
                                    return Arrays.asList(
                                            a.getValue().get(0),
                                            a.getValue().get(1));
                                } else {
                                    return Arrays.asList(
                                            a.getValue().get(1),
                                            a.getValue().get(0));
                                }
                            })
                    .reduce(Arrays.asList(BigDecimal.ZERO, BigDecimal.ZERO),
                            (a, b) -> Arrays.asList(
                                    a.get(0).add(b.get(0)),
                                    a.get(1).add(b.get(1))));
            rotamerCumulativeWeights.put(rotamer,boundedStatWeight);
        }
        return rotamerCumulativeWeights;
    }


    public static void testCumulativeMarginals(Map<String,List<BigDecimal>> residue, String resname, BigDecimal lowerBound, BigDecimal upperBound, Boolean quiet){
        /**
         * Tests to ensure that the marginal distribution over a residue captures the full Boltzmann distribution
         *
         * Something has likely gone wrong if this is not the case (although it is possible that we have not explored nodes for some residues
         * Also I THINK that we have rounding errors because of the necessity of printing things out to file
         *
         * Returns the sum of the marginal partition function for a residue
         */
        List<BigDecimal> totalStatWeight = residue.values().stream()
                .reduce(Arrays.asList(BigDecimal.ZERO, BigDecimal.ZERO),
                        (a, b) -> Arrays.asList(
                                a.get(0).add(b.get(0)),
                                a.get(0).add(b.get(0))));
        if (!quiet) {
            System.out.println(String.format("Residue %s has marginal bounds of [%.8e, %.8e](L,U). Full pfunc -> [%.8e, %.8e](L,U)",
                    resname,
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
            System.err.println(String.format("WARNING: Residue %s marginal distributions have residuals of [%.8f, %.8f](L,U) %% of the full pfunc",
                    resname,
                    percentDiff.get(0)*100,
                    percentDiff.get(1)*100
                    ));
            //TODO: If the marginal value is significantly less than the full partition function, split weight equally among rotamers
            // This is because if we haven't expanded a node, it won't be in the tree, but this means that the statistical weight for
            // That residue's rotamers won't be accounted for - since we don't have info we must split bounds equally... I think
        }
    }

    public static void testCumulativeMarginals(Map<String,List<BigDecimal>> residue, String resname, BigDecimal lowerBound, BigDecimal upperBound){
        testCumulativeMarginals(residue, resname, upperBound, lowerBound, true);
    }
}
