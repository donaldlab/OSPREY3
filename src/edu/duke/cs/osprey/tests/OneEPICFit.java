/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.TermECalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;

/**
 *
 * Fit just one EPIC term
 * This is to improve EPIC by dealing with problem terms
 * 
 * @author mhall44
 */
public class OneEPICFit {
    
    
    
    public static void main(String[] args){
        //Specify desired system and term here
        args = new String[] {"System.cfg","DEE.cfg"};
        int pos1=5;
        int rc1=7;
        int pos2=2;
        int rc2=12;
        //int pos1=3;
        //int rc1=61;
        //int pos2=12;
        //int rc2=76;
        //ROT:3 1 0 12 12 10
        
        RCTuple RCs = new RCTuple(pos1,rc1,pos2,rc2);
        
        
        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(args);//args are configuration files
        cfp.loadData();
        SearchProblem sp = cfp.getSearchProblem();
        
        TermECalculator tec = new TermECalculator(sp.confSpace, sp.shellResidues, 
                true, false, null, new EPICSettings(), false, pos1, pos2);
            
        tec.calcTupleEnergy(RCs);//this might crash on trying to store the fit,
        //but that's OK, just want to see fitting output
    }
    
}
