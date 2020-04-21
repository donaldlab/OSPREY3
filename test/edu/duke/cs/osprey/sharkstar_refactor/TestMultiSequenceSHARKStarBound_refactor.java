package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound;
import edu.duke.cs.osprey.sharkstar.TestSHARKStarBound;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.jetbrains.annotations.NotNull;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

public class TestMultiSequenceSHARKStarBound_refactor extends TestBase {
    private static Molecule metallochaperone;
    private static Molecule protein_1a0r;
    private static PartitionFunctionFactory MSSHARKStarPfuncFactory;

    @BeforeClass
    public static void beforeClass() {
        metallochaperone = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");
        //protein_1a0r = PDBIO.readFile("examples//1a0r_prepped.pdb");
    }

    private MultiSequenceSHARKStarBound_refactor makeSHARKStarPfuncForConfSpace(SimpleConfSpace confSpace, @NotNull Sequence sequence, double epsilon, MultiSequenceSHARKStarBound_refactor preComputedFlex, UpdatingEnergyMatrix preCompCorrections){
        // Set up partition function requirements
        Parallelism parallelism = Parallelism.makeCpu(1);
        ForcefieldParams ffparams = new ForcefieldParams();

        // how should we compute energies of molecules?
        EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();
        // how should we define energies of conformations?
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                        .build()
                        .calcReferenceEnergies()
                )
                .build();

        // BBK* needs rigid energies too
        EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                .setIsMinimizing(false)
                .build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, ecalcRigid);

        EnergyMatrix minEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .setCacheFile(new File("test"+"."+"minimizing"+".emat"))
                .build()
                .calcEnergyMatrix();

        EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .setCacheFile(new File("test"+"."+"minimizing"+".emat"))
                .build()
                .calcEnergyMatrix();

        RCs rcs = sequence.makeRCs(confSpace);

        return new MultiSequenceSHARKStarBound_refactor(confSpace, rigidEmat, minEmat, confEcalcMinimized, rcs, parallelism );

    }

    /**
     * Creates a flexible only confspace with one chain using a few residues from 1CC8
     */
    public SimpleConfSpace make1CC8Flexible(){
        Strand strand1 = new Strand.Builder(metallochaperone).setResidues("A2", "A10").build();
        strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
        strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType);
        strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
        strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType);
        strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType);

        return new SimpleConfSpace.Builder()
                .addStrands(strand1)
                .setShellDistance(9)
                .build();

    }

    @Test
    public void testMakePfunc(){
        SimpleConfSpace confSpace = make1CC8Flexible();
        Sequence wildType = confSpace.makeWildTypeSequence();
        makeSHARKStarPfuncForConfSpace(confSpace, wildType, 0.99, null, null);
    }
}
