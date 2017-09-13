/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;

/**
 *
 * @author Anna Lowegard(anna.lowegard@duke.edu)
 * Adegoke Ojewole (ao68@duke.edu)
 */
public class EWAKRatios {

    public EWAKRatios(ConfigFileParser cfp) {
        // parse config file
    	EWAKConfigFileParser ecfp = new EWAKConfigFileParser(cfp);
        
    	// make search problem
        ecfp.getSearchProblems();
        ecfp.loadEnergyMatrices();
        ecfp.pruneMatrices();
   
    }
    
    public void run() {
    };
    
}
