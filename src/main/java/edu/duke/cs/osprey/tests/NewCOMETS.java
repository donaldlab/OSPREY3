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

import edu.duke.cs.osprey.astar.comets.LME;
import edu.duke.cs.osprey.astar.comets.NewCOMETSDoer;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.gmec.PruningSettings;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author mhall44
 */
public class NewCOMETS {
    
    public static void main(String[] args){
        int numStates = 4;
        SimpleConfSpace[] confSpaces = new SimpleConfSpace[numStates];
        PrecomputedMatrices[] precompMats = new PrecomputedMatrices[numStates];
        LME objFcn = new LME("1 -1 0 0 0",4);
        LME[] constraints = new LME[] {new LME("1 -1 0 0 20",4),new LME("0 0 -1 1 70",4),new LME("0 1 0 1 20",4)};
        int boundMutPos[] = new int[] {3,4,5,6,7};
        int unboundMutPos[] = new int[] {0,1,2,3,4};
        ArrayList<ArrayList<Integer>> mutable2StatePosNums = toDoubleList(boundMutPos,unboundMutPos,boundMutPos,unboundMutPos);
        
        ArrayList<ArrayList<String>> AATypeOptions = toDoubleList(
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ASP","GLU"},
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ASN","GLN","SER","THR"}
        );

        int numMaxMut = -1;
        String wtSeq[] = null;
        int numSeqsWanted = 5;
        ConfEnergyCalculator[] confECalc = new ConfEnergyCalculator[numStates];        
        String stateNames[] = new String[] {"3K75.b","3K75.ub","3LQC.b","3LQC.ub"};
        
        PruningSettings pruningSettings = new PruningSettings();
        pruningSettings.typedep = true;

        boolean useERef = true;//this might seem irrelevant to COMETS
        //but if does affect values of many LME's; the shift in LME value due to eref's
        //is only constant wrt sequence if the LME coefficients sum to 1

        
        for(int state=0; state<numStates; state++){
            confSpaces[state] = prepareConfSpace(state,AATypeOptions);
	    ForcefieldParams ffparams = new ForcefieldParams();
            ffparams.solvScale = 0.;
	    EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces[state], ffparams).build();
            
            ConfEnergyCalculator.Builder confEcalcBuilder = new ConfEnergyCalculator.Builder(confSpaces[state], ecalc);
            if(useERef){
                SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(confSpaces[state],ecalc).build();
                confEcalcBuilder.setReferenceEnergies(eref);
            }
            confECalc[state] = confEcalcBuilder.build();
            
            EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confECalc[state])
				.build()
				.calcEnergyMatrix();
            
            double Ival = 0;
            double Ew = 0;
            precompMats[state] = new PrecomputedMatrices(Ival, Ew, stateNames[state], emat, 
                    confSpaces[state], ecalc, confECalc[state], new EPICSettings(), new LUTESettings(),
                    pruningSettings);//rigid design
        }
            
        NewCOMETSDoer cd = new NewCOMETSDoer(confSpaces,precompMats,objFcn,constraints,
            mutable2StatePosNums,AATypeOptions,numMaxMut,wtSeq,numSeqsWanted,true,confECalc);
        ArrayList<String> bestSequences = cd.calcBestSequences();
    }
    
    
    private static SimpleConfSpace prepareConfSpace(int state, ArrayList<ArrayList<String>> AATypeOptions){

        String pdbFile;
        String[] mutResNums;
        
        switch(state){
            case 0:
                pdbFile = "examples/3K75.3LQC/3K75.b.shell.pdb";
                mutResNums = new String[] {"0391","0409","0411","0422","0424"};
                break;
            case 1:
                pdbFile = "examples/3K75.3LQC/3K75.ub.shell.pdb";
                mutResNums = new String[] {"0291","0309","0311","0322","0324"};
                break;
            case 2:
                pdbFile = "examples/3K75.3LQC/3LQC.b.shell.pdb";
                mutResNums = new String[] {"0591","0609","0611","0622","0624"};
                break;
            case 3:
                pdbFile = "examples/3K75.3LQC/3LQC.ub.shell.pdb";
                mutResNums = new String[] {"0291","0309","0311","0322","0324"};
                break;
            default:
                throw new RuntimeException("Unrecognized state");
        }
        
        Molecule mol = PDBIO.readFile(pdbFile);

        Strand strand = new Strand.Builder(mol).build();
        for(int mutPos=0; mutPos<AATypeOptions.size(); mutPos++)
            strand.flexibility.get(mutResNums[mutPos]).setLibraryRotamers(AATypeOptions.get(mutPos));
        
        if(state%2==0){//bound state, set flexibility for non-designed chain
            strand.flexibility.get("067").setLibraryRotamers("Phe");
            strand.flexibility.get("090").setLibraryRotamers("Thr");
            strand.flexibility.get("0136").setLibraryRotamers("Tyr");
        }
        
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
        return confSpace;
    }
    
    private static ArrayList<ArrayList<Integer>> toDoubleList(int[]... arr){
        ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
        for(int[] a : arr){
            ArrayList<Integer> subAns = new ArrayList<>();
            for(int b : a)
                subAns.add(b);
            ans.add(subAns);
        }
        return ans;
    }
    
    private static ArrayList<ArrayList<String>> toDoubleList(String[]... arr){
        ArrayList<ArrayList<String>> ans = new ArrayList<>();
        for(String[] a : arr){
            ans.add( new ArrayList(Arrays.asList(a)) );
        }
        return ans;
    }
    
}
