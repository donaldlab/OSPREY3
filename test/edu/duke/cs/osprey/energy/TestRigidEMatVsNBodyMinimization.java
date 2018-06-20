package edu.duke.cs.osprey.energy;


import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import java.io.File;

public class TestRigidEMatVsNBodyMinimization {

    public static class ConfSpaces {
        public ForcefieldParams ffparams;
        public SimpleConfSpace protein;
        public SimpleConfSpace ligand;
        public SimpleConfSpace complex;
    }


    @Test
    public void sanityTest() {
        ConfSpaces confSpaces = make1GUASmall(8);
        int[] conf = new int[]{7, 9, 5, 7, 9, 7, 1,4,11};
        Parallelism parallelism = Parallelism.makeCpu(1);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .setIsMinimizing(true)
                .build();
        // Define the rigid energy calculator
        EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
                .setIsMinimizing(false)
                .build();
        // how should we define energies of conformations?
        KStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, minimizingEcalc)
                        .setCacheFile(new File("test.eref.emat"))
                        .build()
                        .calcReferenceEnergies()
                )
                .build();
        SimplerEnergyMatrixCalculator.Builder rigidBuilder = new SimplerEnergyMatrixCalculator.Builder(confEcalcFactory.make(confSpaces.complex, rigidEcalc));
        ConfEnergyCalculator minConfECalc = confEcalcFactory.make(confSpaces.complex, minimizingEcalc);
        SimplerEnergyMatrixCalculator.Builder minimizingBuilder = new SimplerEnergyMatrixCalculator.Builder(minConfECalc);
        EnergyMatrix rigidEmat = rigidBuilder.build().calcEnergyMatrix();
        EnergyMatrix pairwiseMinEmat = minimizingBuilder.build().calcEnergyMatrix();
        PairwiseGScorer pairwiseMinScorer = new PairwiseGScorer(pairwiseMinEmat);
        PairwiseGScorer rigidScorer = new PairwiseGScorer(rigidEmat);
        double rigidScore = rigidScorer.calc(conf);
        double minimizingScore = pairwiseMinScorer.calc(conf);
        ConfSearch.ScoredConf sconf = new ConfSearch.ScoredConf(conf, minimizingScore);
        ConfSearch.EnergiedConf econf = minConfECalc.calcEnergy(sconf);
        System.out.println("Rigid: " + rigidScore + ", pairwise: " + minimizingScore + ", minimized: " + econf.getEnergy());
    }


    public static ConfSpaces make1GUASmall(int numFlex) {

        ConfSpaces confSpaces = new ConfSpaces();

        // configure the forcefield
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

        // make sure all strands share the same template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
                .addMoleculeForWildTypeRotamers(mol)
                .build();

        // define the protein strand
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("1", "180")
                .build();
        int start = 21;
        for (int i = start; i < start + numFlex; i++) {
            protein.flexibility.get(i + "").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

        // define the ligand strand
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("181", "215")
                .build();
        ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers();

        // make the complex conf space ("complex" SimpleConfSpace, har har!)
        confSpaces.protein = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();
        confSpaces.ligand = new SimpleConfSpace.Builder()
                .addStrand(ligand)
                .build();
        confSpaces.complex = new SimpleConfSpace.Builder()
                .addStrands(protein, ligand)
                .build();

        return confSpaces;
    }


}
