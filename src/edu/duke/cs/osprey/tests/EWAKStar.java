/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.ewakstar.NewEWAKStarDoer;
import edu.duke.cs.osprey.confspace.Sequence;
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
import edu.duke.cs.osprey.kstar.KStar;
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
        int numStates = 3;
        Integer[] pos = new Integer[]{3,4,5,6,7};
        SimpleConfSpace[] confSpaces = new SimpleConfSpace[numStates];
        PrecomputedMatrices[] precompMats = new PrecomputedMatrices[numStates];
        ArrayList<Integer> boundMutPos = new ArrayList<> (Arrays.asList(pos));
        
        ArrayList<ArrayList<String>> AATypeOptions = toDoubleList(
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ASP","GLU"},
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ASN","GLN","SER","THR"}
        );

        String wtSeq[] = null;
        int numSeqsWanted = 1;
        ConfEnergyCalculator[] confECalc = new ConfEnergyCalculator[numStates];
        String[] stateNames = new String[] {"3K75.pl","3K75.l","3K75.p"};
        
        PruningSettings pruningSettings = new PruningSettings();
        pruningSettings.typedep = true;
        boolean useERef = true;

        double unboundEw = 1;

        for(int state=0; state<numStates; state++){
            confSpaces[state] = prepareConfSpace(state,AATypeOptions);
            ForcefieldParams ffparams = new ForcefieldParams();
            EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces[state], ffparams).build();

            ConfEnergyCalculator.Builder confEcalcBuilder = new ConfEnergyCalculator.Builder(confSpaces[state], ecalc);
            if(useERef){
                SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(confSpaces[state],ecalc).build();
                confEcalcBuilder.setReferenceEnergies(eref);
            }
            confECalc[state] = confEcalcBuilder.build();

            EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confECalc[state])
                    .setCacheFile(new File("ewak."+Integer.toString(state)+".emat"))
                    .build()
                    .calcEnergyMatrix();

            double Ival = 0;
            double Ew = 0;
            if (state==1) {
                Ew = unboundEw;
            }
            precompMats[state] = new PrecomputedMatrices(Ival, Ew, stateNames[state], emat,
                    confSpaces[state], ecalc, confECalc[state], new EPICSettings(), new LUTESettings(),
                    pruningSettings);//rigid design
        }


        NewEWAKStarDoer ed = new NewEWAKStarDoer(confSpaces,precompMats,
            boundMutPos,AATypeOptions,wtSeq,numSeqsWanted, confECalc, unboundEw);

        ArrayList<String> bestSequences = ed.calcBestSequences();

//        for (int i=0; i<bestSequences.size(); i++) {
//            System.out.println(bestSequences.get(i));
//        }


    }

    private static String printSeqs(ArrayList<String> seqs){
        StringBuffer buf = new StringBuffer();

        for (int i=0; i<seqs.size();i++){
            buf.append(seqs.get(i));
            buf.append('\n');
        }

        return buf.toString();

    }
    
    private static SimpleConfSpace prepareConfSpace(int state, ArrayList<ArrayList<String>> AATypeOptions){

        String pdbFile = "examples/3K75.3LQC/3K75.b.shell.pdb";
        String[] mutResNums;
        String startRes;
        String endRes;

        switch(state){
            case 0:
                startRes = "040";
                endRes = "0428";
                mutResNums = new String[] {"0391","0409","0411","0422","0424"};
                break;
            case 1:
                startRes = "0363";
                endRes = "0428";
                mutResNums = new String[] {"0391","0409","0411","0422","0424"};
                break;
            case 2:
                startRes = "040";
                endRes = "0144";
                mutResNums = new String[] {"067","090","0136"};
                break;
            default:
                throw new RuntimeException("Unrecognized state");
        }

        Molecule mol = PDBIO.readFile(pdbFile);
        Strand strand = new Strand.Builder(mol).setResidues(startRes, endRes).build();

        if(state==0 || state==1){
            for(int mutPos=0; mutPos<AATypeOptions.size(); mutPos++)
                strand.flexibility.get(mutResNums[mutPos]).setLibraryRotamers(AATypeOptions.get(mutPos)).setContinuous();
        }

        if(state==0 || state==2){//bound state, set flexibility for non-designed chain, unbound state, set flexibility
            strand.flexibility.get("067").setLibraryRotamers("Phe").setContinuous();
            strand.flexibility.get("090").setLibraryRotamers("Thr").setContinuous();
            strand.flexibility.get("0136").setLibraryRotamers("Tyr").setContinuous();
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
