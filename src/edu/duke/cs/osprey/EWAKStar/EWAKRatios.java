/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.EWAKStar;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;

/**
 *
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Anna Lowegard(anna.lowegard@duke.edu)
 */
public class EWAKRatios {

    public EWAKRatios(ConfigFileParser cfp) {
        EWAKConfigFileParser ecfp = new EWAKConfigFileParser(cfp);
        SearchProblem[] sps = ecfp.makeSearchProblems();
    }
    
    public void run() {
    };
    
}
