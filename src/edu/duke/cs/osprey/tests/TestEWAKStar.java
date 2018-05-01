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
public class TestEWAKStar {
    
    public static void main(String[] args){
        Integer[] pos = new Integer[]{0,1,2,3,4,5,6,7};
        Integer[] posL = new Integer[]{0,1,2,3,4};
        Integer[] posP = new Integer[]{0,1,2};

        ArrayList<Integer> boundMutPos = new ArrayList<> (Arrays.asList(pos));
        ArrayList<Integer> LMutPos = new ArrayList<> (Arrays.asList(posL));
        ArrayList<Integer> PMutPos = new ArrayList<> (Arrays.asList(posP));

        ArrayList<ArrayList<String>> AATypeOptions = toDoubleList(
                new String[] {"PHE"},
                new String[] {"THR"},
                new String[] {"TYR"},
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
        double unboundEw = 15.0;
        double boundEw = 15.0;
        double ewakstarEw = 1.0;
        double Ival = 0.0;
        String startResPL = "040";
        String endResPL = "0428";
        String startResL = "0363";
        String endResL = "0428";
        String startResP = "040" ;
        String endResP = "0363";
        String pdbFile = "examples/3K75.3LQC/3K75.b.shell.pdb";
        String[] resNumsPL = new String[] {"067","090","0136","0391","0409","0411","0422","0424"};
        String[] resNumsL = new String[]{"0391","0409","0411","0422","0424"};
        String[] resNumsP = new String[]{"067","090","0136"};

        Molecule mol = PDBIO.readFile(pdbFile);
        Strand strandPL = new Strand.Builder(mol).setResidues(startResPL, endResPL).build();
        Strand strandL = new Strand.Builder(mol).setResidues(startResL, endResL).build();
        Strand strandP = new Strand.Builder(mol).setResidues(startResP, endResP).build();

        for(Integer p: pos) {
            strandPL.flexibility.get(resNumsPL[p]).setLibraryRotamers(AATypeOptions.get(p)).setContinuous();
        }
        for(Integer p: posL) {
            strandL.flexibility.get(resNumsL[p]).setLibraryRotamers(AATypeOptions.get(p)).setContinuous();
        }
        for(Integer p: posP) {
            strandP.flexibility.get(resNumsP[p]).setLibraryRotamers(AATypeOptions.get(p)).setContinuous();
        }

        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strandPL).build();
        SimpleConfSpace confSpaceP = new SimpleConfSpace.Builder().addStrand(strandP).build();
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

        PrecomputedMatrices precompMat = new PrecomputedMatrices(Ival, boundEw, "PL", emat,
                confSpace, ecalc, confECalc, new EPICSettings(), new LUTESettings(),
                pruningSettings);

        String LmatrixName = "ewak.L.emat";

        String PmatrixName = "ewak.P.emat";

        NewEWAKStarDoer ed = new NewEWAKStarDoer(confSpace, confSpaceL, confSpaceP, precompMat,
            boundMutPos, LMutPos, PMutPos, AATypeOptions, numSeqsWanted, confECalc, orderOfMag, unboundEw, boundEw,
                ewakstarEw, startResL, endResL, startResP, endResP, mol, resNumsPL, resNumsL, resNumsP, Ival, pruningSettings,
                LmatrixName, PmatrixName, ffparams);


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
