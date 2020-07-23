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
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
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
        System.out.println(String.format("Problematic triple: %s, LB: %.9f, Min E: %.9f",
                problematicTriple.toString(),
                emat.getInternalEnergy(problematicTriple),
                confEcalc.calcEnergy(problematicTriple).energy
                ));
        for (ConfSearch.ScoredConf conf : confs){
            System.out.println(conf);
            System.out.println(confEcalc.calcEnergy(conf));

        }
    }
}
