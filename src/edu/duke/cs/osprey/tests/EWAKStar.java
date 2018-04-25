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
        int numStates = 2;
        Integer[] pos = new Integer[]{3,4,5,6,7};
        Integer[] posL = new Integer[]{0,1,2,3,4};
        ArrayList<Integer> boundMutPos = new ArrayList<> (Arrays.asList(pos));
        ArrayList<Integer> unboundMutPos = new ArrayList<> (Arrays.asList(posL));
        
        ArrayList<ArrayList<String>> AATypeOptions = toDoubleList(
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ASP","GLU"},
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ILE","LEU","MET","PHE","TRP","TYR","VAL"},
                new String[] {"ASN","GLN","SER","THR"}
        );

        int numSeqsWanted = 10000;
        
        PruningSettings pruningSettings = new PruningSettings();
        pruningSettings.typedep = true;

        double orderOfMag = 5.0;
        double unboundPFw = 15.0;
        double boundPFw = 15.0;
        double Ival = 0.0;
        String startResPL = "040";
        String endResPL = "0428";
        String startResL = "0363";
        String endResL = "0428";
        String pdbFile = "examples/3K75.3LQC/3K75.b.shell.pdb";
        String[] mutResNums = new String[] {"0391","0409","0411","0422","0424"};

        Molecule mol = PDBIO.readFile(pdbFile);
        Strand strandPL = new Strand.Builder(mol).setResidues(startResPL, endResPL).build();
        Strand strandL = new Strand.Builder(mol).setResidues(startResL, endResL).build();

        for(int mutPos=0; mutPos<AATypeOptions.size(); mutPos++) {
            strandPL.flexibility.get(mutResNums[mutPos]).setLibraryRotamers(AATypeOptions.get(mutPos)).setContinuous();
            strandL.flexibility.get(mutResNums[mutPos]).setLibraryRotamers(AATypeOptions.get(mutPos)).setContinuous();
        }

        strandPL.flexibility.get("067").setLibraryRotamers("Phe").setContinuous();
        strandPL.flexibility.get("090").setLibraryRotamers("Thr").setContinuous();
        strandPL.flexibility.get("0136").setLibraryRotamers("Tyr").setContinuous();

        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strandPL).build();
        SimpleConfSpace confSpaceL = new SimpleConfSpace.Builder().addStrand(strandL).build();

        ForcefieldParams ffparams = new ForcefieldParams();
        EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build();

        ConfEnergyCalculator.Builder confEcalcBuilder = new ConfEnergyCalculator.Builder(confSpace, ecalc);

        //use reference energies
        SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(confSpace,ecalc).build();
        confEcalcBuilder.setReferenceEnergies(eref);

        ConfEnergyCalculator confECalc = confEcalcBuilder.build();

        EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confECalc)
                .setCacheFile(new File("ewak.PL.emat"))
                .build()
                .calcEnergyMatrix();

        PrecomputedMatrices precompMat = new PrecomputedMatrices(Ival, boundPFw, "PL", emat,
                confSpace, ecalc, confECalc, new EPICSettings(), new LUTESettings(),
                pruningSettings);

        String LmatrixName = "ewak.L.emat";

        NewEWAKStarDoer ed = new NewEWAKStarDoer(confSpace, confSpaceL, precompMat,
            boundMutPos, unboundMutPos, AATypeOptions, numSeqsWanted, confECalc, orderOfMag, unboundPFw, boundPFw,
                startResL, endResL, mol, mutResNums, Ival, pruningSettings, LmatrixName);

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

    
    private static ArrayList<ArrayList<String>> toDoubleList(String[]... arr){
        ArrayList<ArrayList<String>> ans = new ArrayList<>();
        for(String[] a : arr){
            ans.add( new ArrayList(Arrays.asList(a)) );
        }
        return ans;
    }
    
}
