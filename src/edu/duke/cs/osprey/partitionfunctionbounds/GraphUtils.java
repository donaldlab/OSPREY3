/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import Jama.*;
import java.math.BigDecimal;
import java.math.RoundingMode;

/**
 *
 * @author hmn5
 */
public class GraphUtils {
	
    boolean[][] interactionGraph;

    public GraphUtils(boolean[][] interactionGraph) {
        this.interactionGraph = interactionGraph;
    }
        
    private static Matrix getAdjacencyMatrix(boolean[][] interactionGraph) {
        double[][] adjMatAsList = new double[interactionGraph.length][];
        for (int i = 0; i < interactionGraph.length; i++) {
            adjMatAsList[i] = new double[interactionGraph.length];
            for (int j = 0; j < i; j++) {
                adjMatAsList[i][j] = interactionGraph[i][j] ? 1. : 0.;
                adjMatAsList[j][i] = adjMatAsList[i][j];
            }
        }
        Matrix adjMat = new Matrix(adjMatAsList);
        return adjMat;
    }

    /**
     *
     * @param edgeWeights Just contains the upper-triangular information of the
     * symmetric matrix. It must be converted into a symmetric matrix
     * @return
     */
    private static Matrix getWeightedAdjMatrix(double[][] edgeWeights) {
        double[][] weightedAdjMatAsList = new double[edgeWeights.length][edgeWeights.length];
        for (int i = 0; i < edgeWeights.length; i++) {
            for (int j = 0; j < i; j++) {
                weightedAdjMatAsList[i][j] = Math.abs(edgeWeights[i][j]);
                weightedAdjMatAsList[j][i] = weightedAdjMatAsList[i][j];
            }
        }
        Matrix weightAdjMat = new Matrix(weightedAdjMatAsList);
        return weightAdjMat;
    }

    private static Matrix getLaplacianMatrix(Matrix adjacentMatrix) {
        Matrix laplacian = adjacentMatrix.copy().uminus();
        for (int i = 0; i < adjacentMatrix.getColumnDimension(); i++) {
            double degreeOfVertexI = 0;
            for (int j = 0; j < adjacentMatrix.getRowDimension(); j++) {
                degreeOfVertexI += adjacentMatrix.get(i, j);
            }
            laplacian.set(i, i, degreeOfVertexI);
        }
        return laplacian;
    }

    private static Matrix subtractOneElementwise(Matrix matrix) {
        double[][] ones = new double[matrix.getColumnDimension()][];
        for (int i = 0; i < matrix.getColumnDimension(); i++) {
            ones[i] = new double[matrix.getRowDimension()];
            for (int j = 0; j < matrix.getRowDimension(); j++) {
                ones[i][j] = 1d;
            }
        }
        Matrix onesMatrix = new Matrix(ones);
        return matrix.minus(onesMatrix);
    }

    public static double[][] getEdgeProbabilities(boolean[][] interactionGraph) {
        Matrix adj = getAdjacencyMatrix(interactionGraph);
        Matrix laplacian = getLaplacianMatrix(adj);
        Matrix laplMinusOne = subtractOneElementwise(laplacian);

        Matrix inverseLaplacian = laplMinusOne.inverse();

        double[][] edgeProbs = new double[interactionGraph.length][];
        for (int i = 0; i < interactionGraph.length; i++) {
            edgeProbs[i] = new double[i];
            for (int j = 0; j < i; j++) {
                double prob = adj.get(i, j) * (inverseLaplacian.get(i, i) + inverseLaplacian.get(j, j) - 2 * inverseLaplacian.get(i, j));
                edgeProbs[i][j] = prob;
            }
        }

        return edgeProbs;
    }

    /**
     *
     * @param edgeWeights Just contains the upper-triangular information of the
     * symmetric matrix. It must be converted into a symmetric matrix
     * @return
     */
    public static double[][] getEdgeProbabilities(double[][] edgeWeights) {
        Matrix adj = getWeightedAdjMatrix(edgeWeights);
        Matrix laplacian = getLaplacianMatrix(adj);
        Matrix laplMinusOne = subtractOneElementwise(laplacian);

        Matrix inverseLaplacian = laplMinusOne.inverse();

        double[][] edgeProbs = new double[edgeWeights.length][];
        for (int i = 0; i < edgeWeights.length; i++) {
            edgeProbs[i] = new double[i];
            for (int j = 0; j < i; j++) {
                double prob = adj.get(i, j) * (inverseLaplacian.get(i, i) + inverseLaplacian.get(j, j) - 2 * inverseLaplacian.get(i, j));
                edgeProbs[i][j] = prob;
            }
        }

        return edgeProbs;
    }

    public static double[] getWeightedDegrees(double[][] edgeWeights) {
        Matrix adj = getWeightedAdjMatrix(edgeWeights);
        double[] degrees = new double[edgeWeights.length];
        for (int i = 0; i < adj.getColumnDimension(); i++) {
            for (int j = 0; j < adj.getRowDimension(); j++) {
                degrees[i] += adj.get(i,j);
            }
        }
        return degrees;
    }

    public static double round(double value, int places) {
        if (places < 0) {
            throw new IllegalArgumentException();
        }

        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

    /*    private Matrix getWeightedAdjMatrix(boolean[][] interactionGraph,
     ArrayList<MRFNode> nodeList, UpdatedEmat emat) {
     double[][] adjMatAsList = new double[nodeList.size()][];
     for (int i = 0; i < nodeList.size(); i++) {
     MRFNode nodeI = nodeList.get(i);
     adjMatAsList[i] = new double[nodeList.size()];
     for (int j = 0; j < i; j++) {
     MRFNode nodeJ = nodeList.get(j);
     if (!interactionGraph[i][j]) {
     adjMatAsList[i][j] = 0.;
     adjMatAsList[j][i] = adjMatAsList[i][j];
     } else {
     adjMatAsList[i][j] = getEdgeWeight(nodeI, nodeJ, emat);
     }
     }
     }
     return 
     }
     */
}
