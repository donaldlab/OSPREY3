package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;

import static edu.duke.cs.osprey.sharkstar.TestSHARKStar.loadFromCFS;

public class TestConfEnergyCalculator extends TestBase {

    private static SimpleConfSpace confSpace_4wyq_E_complex;
    private static SimpleConfSpace confSpace_4wyq_E_protein;
    private static ForcefieldParams ffparams;
    private static List<ConfSearch.ScoredConf> confs;
    private static EnergyMatrix emat;
    private static ConfEnergyCalculator confEcalc;

    private static RCTuple problematicTriple = new RCTuple(0,15,2,5,4,12);

    @BeforeClass
    public static void beforeClass() {
        // get a conf space
        Molecule mol = PDBIO.readFile("examples/python.KStar/4wyq_prepped.pdb");
        Strand strand1 = new Strand.Builder(mol)
                .setResidues("D269", "D391")
                .build();
        Strand strand2 = new Strand.Builder(mol)
                .setResidues("E289", "E363")
                .build();
        strand1.flexibility.get("D285").setLibraryRotamers("ASP").addWildTypeRotamers().setContinuous();
        strand2.flexibility.get("E340").setLibraryRotamers("SER").addWildTypeRotamers().setContinuous();
        strand2.flexibility.get("E324").setLibraryRotamers("GLN").addWildTypeRotamers().setContinuous();
        strand2.flexibility.get("E318").setLibraryRotamers("ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIP", "HIE", "HID", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
        strand2.flexibility.get("E322").setLibraryRotamers("LEU").addWildTypeRotamers().setContinuous();
        strand2.flexibility.get("E320").setLibraryRotamers("ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIP", "HIE", "HID", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
        //confSpace_4wyq_E_complex = new SimpleConfSpace.Builder().addStrands(strand1, strand2).build();
        confSpace_4wyq_E_protein = new SimpleConfSpace.Builder().addStrand(strand2).build();

        ffparams = new ForcefieldParams();

        // get an energy matrix
        EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace_4wyq_E_protein, ffparams)
                .setParallelism(Parallelism.makeCpu(4))
                .setIsMinimizing(true)
                .build();
        confEcalc= new ConfEnergyCalculator.Builder(confSpace_4wyq_E_protein, ecalc)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace_4wyq_E_protein, ecalc)
                        .build()
                        .calcReferenceEnergies()
                )
                .build();
        emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
        .build()
        .calcEnergyMatrix();

        Function<int[], ConfSearch.ScoredConf> maker = (conf) -> new ConfSearch.ScoredConf(conf, emat.confE(conf));

        confs = new ArrayList<>();
        confs.add(maker.apply(new int[] {15,198,5,4,12}));
        confs.add(maker.apply(new int[] {15,35,5,4,12}));
        confs.add(maker.apply(new int[] {15,201,5,4,12}));
    }

    @Test
    public void testLowerBounds(){
        ConfAnalyzer analyzer = new ConfAnalyzer(confEcalc);
        double tripleLB = emat.getInternalEnergy(problematicTriple);
        double tripleE = confEcalc.calcEnergy(problematicTriple).energy;
        double tripleCorrection = tripleE-tripleLB;
        System.out.println(String.format("Problematic triple: %s, LB: %.9f, Min E: %.9f, Correction %.9f",
                problematicTriple.toString(),
                tripleLB,
                tripleE,
                tripleCorrection
                ));
        for (ConfSearch.ScoredConf conf : confs){
            System.out.println(conf);
            System.out.println(confEcalc.calcEnergy(conf));
            ConfAnalyzer.ConfAnalysis analysis = analyzer.analyze(conf);
            EnergyMatrix energyAnalysis = analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All);
            EnergyMatrix scoreAnalysis = analysis.breakdownScoreByPosition(emat);
            EnergyMatrix diff = energyAnalysis.diff(scoreAnalysis);

            /*
            for (int pos = 0; pos < diff.getNumPos(); pos++) {
                for (int rc = 0; rc < diff.getNumConfAtPos(pos); rc++) {
                    for (int pos2 = 0; pos2 < diff.getNumPos(); pos2++) {
                        for (int rc2 = 0; rc2 < diff.getNumConfAtPos(pos2); rc2++) {
                            if (pos >= pos2)
                                continue;
                            double sum = 0;
                            sum += diff.getOneBody(pos, rc);
                            sum += diff.getPairwise(pos, rc, pos2, rc2);
                            sum += diff.getOneBody(pos2, rc2);
                            System.out.println(sum);
                        }
                    }
                }
            }
             */
            double diffSum = 0;
            for(int i : problematicTriple.pos){
                diffSum += diff.getOneBody(i,0);
                for (int j : problematicTriple.pos){
                    if (i>=j)
                        continue;
                    diffSum += diff.getPairwise(i,0, j, 0);
                }
            }
            System.out.println(String.format("Diff Sum: %f",diffSum));
            System.out.println(String.format("Corrected LB: %f", conf.getScore()+ tripleCorrection));

        }
    }

    @Test
    /**
     * These are conformations for which I am seeing pairwise-minimized lower bounds that do not properly bound minimized energy.
     * The behavior looks non-deterministic, as repeated running of this test results in different outcomes.
     */
    public void problemConfs4wyuA(){
        // get a conf space
        Molecule mol = PDBIO.readFile("examples/python.KStar/4wyu_prepped.pdb");
        Strand strand1 = new Strand.Builder(mol)
                .setResidues("D-6", "D0")
                .build();
        Strand strand0 = new Strand.Builder(mol)
                .setResidues("A3", "A203")
                .build();
        strand1.flexibility.get("D-4").setLibraryRotamers(Strand.WildType, "PHE").addWildTypeRotamers().setContinuous();
        strand0.flexibility.get("A176").setLibraryRotamers(Strand.WildType, "LEU").addWildTypeRotamers().setContinuous();
        strand0.flexibility.get("A177").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIP", "HIE", "HID", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand0).build();

        ForcefieldParams ffparams = new ForcefieldParams();

        /*
        try {
            confSpace = loadFromCFS("test-resources/4wyu_A_3res_6.090E+03.cfs").protein;
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
         */

        // get an energy matrix
        EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(Parallelism.makeCpu(4))
                .setIsMinimizing(true)
                .build();
        ConfEnergyCalculator confEcalc= new ConfEnergyCalculator.Builder(confSpace, ecalc)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
                        .build()
                        .calcReferenceEnergies()
                )
                .build();
        EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
                .build()
                .calcEnergyMatrix();

        ConfAnalyzer analyzer = new ConfAnalyzer(confEcalc);

        List<int[]> confs = new ArrayList<>();
        confs.add(new int[]{1,8});
        confs.add(new int[]{3,6});
        confs.add(new int[]{3,60});
        confs.add(new int[]{3,62});
        confs.add(new int[]{3,152});
        confs.add(new int[]{3,155});
        confs.add(new int[]{5,155});
        confs.add(new int[]{3,168});
        confs.add(new int[]{3,176});
        confs.add(new int[]{4,176});


        for (int[] conf : confs){
            double lb = emat.confE(conf);
            //double energy = analyzer.analyze(new ConfSearch.ScoredConf(conf, lb)).epmol.energy;
            double energy = confEcalc.calcEnergy(new ConfSearch.ScoredConf(conf, lb)).getEnergy();
            System.out.println(String.format("For %s, test whether LB: %f <= E: %f --> %b",
                    Arrays.toString(conf),
                    lb,
                    energy,
                    lb <= energy
                    ));
            //assert(lb < energy);
        }

    }
}
