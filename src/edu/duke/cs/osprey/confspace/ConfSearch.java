/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.math.BigInteger;

//This is a general interface for things that search conformational space
//like an A* trees, a BWM* search, or a WCSP solver
//For each of these, we instantiate the ConfSearch, 
//call nextConf to get the GMEC for the model being searched,
//and then successive calls to nextConf() return the next conformation in the gap-free list
/**
 *
 * @author mhall44
 */
public interface ConfSearch {
    
	BigInteger getNumConformations();
    int[] nextConf();
    
}
