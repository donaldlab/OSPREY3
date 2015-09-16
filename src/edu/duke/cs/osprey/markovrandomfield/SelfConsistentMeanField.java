/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.markovrandomfield;

import java.math.BigDecimal;
import java.util.ArrayList;
/**
 *
 * @author hmn5
 */
public class SelfConsistentMeanField implements InferenceCalculator{
    
    BigDecimal partitionFunction;
    
    ArrayList<MRFNode> nodeList;
    
    final int maxNumberIterations = 10000;
    public SelfConsistentMeanField(MarkovRandomField mrf){
        
    }
    
    
    @Override
    public BigDecimal calcPartitionFunction(){
        BigDecimal partitionFunction = new BigDecimal(1.0);
        //...
        return partitionFunction;
    }
}
