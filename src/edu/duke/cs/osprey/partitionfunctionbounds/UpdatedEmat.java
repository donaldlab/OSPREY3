/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.confspace.TupleMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class UpdatedEmat extends TupleMatrix<Double> {

    //List of clamped nodes (ie. only have one label)
    ArrayList<MRFNode> clampedNodeList;
    boolean[][] interactionGraph;

    //The const term will be the original const term plus the internal energies
    //of all the clamped nodes
    private double consTerm;

    public UpdatedEmat(EnergyMatrix emat, ArrayList<MRFNode> aClampedNodeList, boolean[][] aInteractionGraph) {
        super(emat);
        this.consTerm = emat.getConstTerm();
        this.clampedNodeList = aClampedNodeList;
        for (MRFNode node : clampedNodeList) {
            if (node.getNumLabels() > 1) {
                throw new RuntimeException("Updated Emat ERROR: Clamped Nodes can only have one label");
            }
            this.consTerm += emat.getOneBody(node.posNum, node.labels[0]);
        }
        for (int i = 0; i < clampedNodeList.size(); i++) {
            for (int j = 0; j < i; j++) {
                MRFNode nodeI = clampedNodeList.get(i);
                MRFNode nodeJ = clampedNodeList.get(j);
                this.consTerm += emat.getPairwise(nodeI.posNum, nodeI.labels[0], nodeJ.posNum, nodeJ.labels[0]);
            }
        }
        this.interactionGraph = aInteractionGraph;
    }

    @Override
    public Double getOneBody(int res, int index) {
        double oneBodyE = this.oneBody.get(res).get(index);
        //add the pairwise energy with clamped nodes;
        for (MRFNode node : clampedNodeList) {
            if (interactionGraph[res][node.posNum]) {
                oneBodyE += this.getPairwise(res, index, node.posNum, node.labels[0]);
            }
        }
        return oneBodyE;
    }

    public double getConstTerm() {
        return this.consTerm;
    }

}
