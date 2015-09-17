/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.markovrandomfield;

import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import java.math.BigDecimal;
import java.util.ArrayList;
/**
 *
 * @author hmn5
 */
public class SelfConsistentMeanField implements InferenceCalculator{
    
    BigDecimal partitionFunction;
    
    ArrayList<MRFNode> nodeList;
    
    EnergyMatrix emat;
    
    final int maxNumberIterations = 10000;
    double constRT = PoissonBoltzmannEnergy.constRT;
    public SelfConsistentMeanField(MarkovRandomField mrf){
        this.nodeList = mrf.nodeList;
        this.emat = mrf.emat;
    }
    
    //Initialize the beliefs of each label in a node to be uniform among all possible
    //labels that belong to the node
    private void initializeBeliefs(){
        for (MRFNode node : this.nodeList){
            int numLabels = node.labelList.size();
            for (MRFLabel label : node.labelList){
                label.currentBelief = 1.0/(double)numLabels;
                label.oldBelief = 1.0/(double)numLabels;
            }
        }
    }
    
    
    //Returns the MeanFieldEnergy for a label
    private double getMeanFieldEnergy(MRFNode node,MRFLabel label){
        double meanFieldE = 0.0;
        for (MRFNode neighbor : node.neighborList){
            //get the average energy between neighbor and our label
            double averageE = getAverageEnergy(node, label, neighbor);
            meanFieldE += averageE;
        }
        return meanFieldE;
    }
    
    //Get average energy of label1 from node1, with respect to all labels of node2
    //this is used as a subroutine in getMeanField()
    private double getAverageEnergy(MRFNode node1, MRFLabel label1, MRFNode node2){
        double averageE = 0.0;
        for (MRFLabel label2 : node2.labelList){
            double E = this.emat.getPairwise(node1.nodeNum, label1.labelNum, node2.nodeNum, label2.labelNum);
            averageE += E*label2.currentBelief;
        }
        return averageE;
    }
    
    //update all of the beliefs belonging to node
    private void updateNodeBeliefs(MRFNode node){
        //create a normalizing constant to normalize beliefs
        BigDecimal partFunction = new BigDecimal(0.0);
        for (MRFLabel label : node.labelList){
            double oneBodyE = this.emat.getOneBody(node.nodeNum, label.labelNum);
            double meanFieldE = getMeanFieldEnergy(node, label);
            //unnormalized updateBelief
            BigDecimal updateBelief
        }
        
    }
    
    @Override
    public BigDecimal calcPartitionFunction(){
        BigDecimal partitionFunction = new BigDecimal(1.0);
        //...
        return partitionFunction;
    }
}
