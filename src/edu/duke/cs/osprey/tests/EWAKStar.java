/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.astar.ewakstar.NewEWAKStarDoer;
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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author lowegard, based on NewCOMETS
 */
public class EWAKStar {
    
    public static void main(String[] args){
        Integer[] pos = new Integer[]{3,4,5,6,7};
        ArrayList<Integer> boundMutPos = new ArrayList<> (Arrays.asList(pos));
        
        ArrayList<ArrayList<String>> AATypeOptions = toDoubleList(
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ASP","GLU"},
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ASN","GLN","SER","THR"}
        );

        String wtSeq[] = null;
        int numSeqsWanted = 6;
        String stateName = "3K75.b";
        
        PruningSettings pruningSettings = new PruningSettings();
        pruningSettings.typedep = true;

        SimpleConfSpace confSpace = prepareConfSpace(AATypeOptions);
	    ForcefieldParams ffparams = new ForcefieldParams();
            ffparams.solvScale = 0.95;
	    EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build();
            
	    ConfEnergyCalculator.Builder confEcalcBuilder = new ConfEnergyCalculator.Builder(confSpace, ecalc);

	    SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(confSpace,ecalc).build();
	    confEcalcBuilder.setReferenceEnergies(eref);

	    ConfEnergyCalculator confECalc = confEcalcBuilder.build();
            
	    EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confECalc)
                .setCacheFile(new File("ewak_test.emat"))
                .build()
				.calcEnergyMatrix();
            
	    double Ival = 0;
	    double Ew = 0;
	    PrecomputedMatrices precompMats = new PrecomputedMatrices(Ival, Ew, stateName, emat,
                    confSpace, ecalc, confECalc, new EPICSettings(), new LUTESettings(),
                    pruningSettings);
            
        NewEWAKStarDoer ed = new NewEWAKStarDoer(confSpace,precompMats,
            boundMutPos,AATypeOptions,wtSeq,numSeqsWanted, confECalc);

        ArrayList<String> bestSequences = ed.calcBestSequences();

        System.out.println(printSeqs(bestSequences));

    }

    private static String printSeqs(ArrayList<String> seqs){
        StringBuffer buf = new StringBuffer();

        for (int i=0; i<seqs.size();i++){
            buf.append(seqs.get(i));
            buf.append('\n');
        }

        return buf.toString();

    }
    
    private static SimpleConfSpace prepareConfSpace(ArrayList<ArrayList<String>> AATypeOptions){

        String pdbFile;
        String[] mutResNums;

        pdbFile = "examples/3K75.3LQC/3K75.b.shell.pdb";
        mutResNums = new String[] {"0391","0409","0411","0422","0424"};

        Molecule mol = PDBIO.readFile(pdbFile);

        Strand strand = new Strand.Builder(mol).build();
        for(int mutPos=0; mutPos<AATypeOptions.size(); mutPos++)
            strand.flexibility.get(mutResNums[mutPos]).setLibraryRotamers(AATypeOptions.get(mutPos)).setContinuous();
        
        //bound state, set flexibility for non-designed chain
        strand.flexibility.get("067").setLibraryRotamers("Phe").setContinuous();
        strand.flexibility.get("090").setLibraryRotamers("Thr").setContinuous();
        strand.flexibility.get("0136").setLibraryRotamers("Tyr").setContinuous();

        
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
