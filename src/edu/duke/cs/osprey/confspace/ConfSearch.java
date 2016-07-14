/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.math.BigInteger;
import java.util.List;

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
    ScoredConf nextConf();
    List<ScoredConf> nextConfs(double maxEnergy);
    
    public static class ScoredConf {
        
        private int[] assignments;
        private double score;
        
        public ScoredConf(int[] assignments, double score) {
            this.assignments = assignments;
            this.score = score;
        }
        
        public ScoredConf(ScoredConf other) {
        	this.assignments = other.assignments.clone();
        	this.score = other.score;
        }
        
        public int[] getAssignments() {
            return assignments;
        }
        
        public double getScore() {
            return score;
        }
    }
    
    public static class EnergiedConf extends ScoredConf {
        
        private double energy;
        
        public EnergiedConf(ScoredConf conf, double energy) {
        	super(conf);
            this.energy = energy;
        }
        
        public EnergiedConf(int[] conf, double score, double energy) {
        	super(conf, score);
        	this.energy = energy;
        }
        
        public double getEnergy() {
            return energy;
        }
    }
}
