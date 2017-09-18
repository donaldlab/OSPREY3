/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.GMECFinder;

/**
 *
 * @author Anna Lowegard(anna.lowegard@duke.edu)
 * Adegoke Ojewole (ao68@duke.edu)
 */
public class EWAKRatios {

    public EWAKRatios(ConfigFileParser cfp) {
        
    	GMECFinder complexes = new GMECFinder();
    	complexes.init(cfp);
    	List<EnergiedConf> complexConfs = complexes.calcGMEC();
    	
    	/* TODO:
    	 * 1) iterate through list of complex confs
    	 * 		map confs to sequence
    	 * 		update partition function of sequence
    	 * 2) make list of P and L only sequences from complex sequences
    	 * 		limit P and L pruning matrices accordingly
    	 * 3) repeat step 1 for P and L
    	 * 4) print output to file
    	 */
    	
    		
    	
    	/*
    	// parse config file
    	EWAKConfigFileParser ecfp = new EWAKConfigFileParser(cfp);
    	// make search problem
        SearchProblem[] sps = ecfp.getSearchProblems();
        ecfp.loadEnergyMatrices();
        ecfp.pruneMatrices();
        
        GMECFinder[] gfs = new GMECFinder[sps.length];
        for(int i = 0; i < gfs.length; ++i) {
        	gfs[i] = new GMECFinder();
        	gfs[i].init(cfp, sps[i]);
        	gfs[i].calcGMEC();
        }
    	*/
    }
    
    public void run() {
    };
    
}
