/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import java.math.BigDecimal;
/**
 *
 * @author hmn5
 */
public interface InferenceCalculator {
    //Any method that provides inference over MRFs (i.e., SCFM) will do so by 
    //calculating a partition funciton.
    
    public abstract BigDecimal calcPartitionFunction();
}
