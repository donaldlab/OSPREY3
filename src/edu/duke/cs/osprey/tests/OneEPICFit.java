/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
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
