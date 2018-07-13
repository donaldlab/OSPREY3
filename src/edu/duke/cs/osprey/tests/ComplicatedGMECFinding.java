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

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.CATSStrandFlex;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.DEEPerStrandFlex;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.confspace.StrandFlex;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.DEEGMECFinder;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * Based on test/GMEC/TestSimpleGMECFinder,
 * but with added complications like DEEPer and LUTE
 * 
 * @author mhall44
 */
public class ComplicatedGMECFinding {
    
    static boolean useDEEGMECFinder = true;
    static boolean useERef = false;
    static boolean useLittleDEEPer = false;
    static boolean useLUTE = true;
    static boolean useEPIC = false;
		
        public static void main(String args[]){
            SimpleConfSpace confSpace = useLittleDEEPer ? littleDEEPerConfSpace() : prepareConfSpace();
	    ForcefieldParams ffparams = new ForcefieldParams();
            ffparams.solvScale = 0.;
	    EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build();
            
            
            ConfEnergyCalculator.Builder confEcalcBuilder = new ConfEnergyCalculator.Builder(confSpace, ecalc);
            if(useERef){
                SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(confSpace,ecalc).build();
                confEcalcBuilder.setReferenceEnergies(eref);
            }
            ConfEnergyCalculator confEcalc = confEcalcBuilder.build();
            
            
            EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build()
				.calcEnergyMatrix();
            
            
            if(useDEEGMECFinder){
                DEEGMECFinder gf = new DEEGMECFinder.Builder( emat, confSpace, ecalc, confEcalc, "1CC8.dee" ).build();
                
                if(useEPIC){
                    gf.epicSettings = EPICSettings.defaultEPIC();
                    gf.pruningSettings.algOption = 3;
                }
                if(useLUTE){
                    gf.luteSettings = LUTESettings.defaultLUTE();
                    gf.pruningSettings.algOption = 3;//nice to prune more if we're going to do LUTE
                    gf.pruningSettings.useTriples = true;
                }
                
                gf.calcGMEC();
            }
            else {
                ConfSearch search = new ConfAStarTree.Builder(emat, confSpace).build();//try regular A* tree next

                SimpleGMECFinder gf = new SimpleGMECFinder.Builder( search, confEcalc ).build();
                gf.find();
            }
        }
        
        private static SimpleConfSpace prepareConfSpace(){
            Molecule mol = PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb");
            String allowedAAs[][] = new String[][] {
                new String[] {"Ile","Ala","Gly"},
                new String[] {"Ser","Ala","Gly"},
                new String[] {"Met","Ser","Ala","Gly"},
                new String[] {"Glu","Ala","Gly"},
                new String[] {"Ala","Gly"},
                new String[] {"Gln","Ala","Gly"},
                new String[] {"Leu","Ala","Gly"} };
            
            String flexResNums[] = new String[] {"38","39","40","41","42","43","44"};
            
            
            Strand strand = new Strand.Builder(mol).build();
            for(int pos=0; pos<7; pos++)
                strand.flexibility.get(flexResNums[pos]).setLibraryRotamers(allowedAAs[pos]).setContinuous();
            
            
            /*StrandFlex bbFlex = new DEEPerStrandFlex( strand, new DEEPerSettings(true, "aa.pert", 
                true, "none", false, 2.5, 2.5, false, new ArrayList(Arrays.asList(flexResNums)), 
                "examples/python.GMEC/1CC8.ss.pdb", false, strand.templateLib) );*/
            //StrandFlex bbFlex = new CATSStrandFlex(strand, "38", "44");
            //just give it the flex res for that strand
            
            
            SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand/*, bbFlex*/).build();
            return confSpace;
        }
        
        
        private static SimpleConfSpace littleDEEPerConfSpace(){
            Molecule mol = PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb");
            String allowedAAs[][] = new String[][] {
                new String[] {"Ile","Ala","Gly"},
                new String[] {"Ser","Ala","Gly"},
                new String[] {"Met","Ser","Ala","Gly","Leu"},
                new String[] {"Glu","Ala","Gly"} };
            
            String flexResNums[] = new String[] {"38","39","40","41"};
            
            
            Strand strand = new Strand.Builder(mol).build();
            for(int pos=0; pos<4; pos++)
                strand.flexibility.get(flexResNums[pos]).setLibraryRotamers(allowedAAs[pos]).setContinuous();
            
            
            StrandFlex bbFlex = new DEEPerStrandFlex( strand, new DEEPerSettings(true, "1CC8.d.pert", 
                true, "none", false, 2.5, 2.5, false, new ArrayList(Arrays.asList(flexResNums)), 
                "examples/python.GMEC/1CC8.ss.pdb", false, strand.templateLib) );
            //just give it the flex res for that strand
            
            
            SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand, bbFlex).build();
            return confSpace;
        }
        
}
