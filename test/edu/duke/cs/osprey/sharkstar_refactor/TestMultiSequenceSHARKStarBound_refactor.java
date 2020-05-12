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
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
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

    private GradientDescentPfunc makeGradientDescentPfuncForConfSpace(SimpleConfSpace confSpace, @NotNull Sequence sequence, double epsilon){
        // Set up partition function requirements
        Parallelism parallelism = Parallelism.makeCpu(4);
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

        PartitionFunctionFactory pfuncFactory = new PartitionFunctionFactory(confSpace, confEcalcMinimized, "pfunc");
        // filter the global sequence to this conf space
        // make the partition function
        RCs rcs = sequence.makeRCs(confSpace);

        return (GradientDescentPfunc) pfuncFactory.makePartitionFunctionFor(rcs, rcs.getNumConformations(), epsilon);
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
                .setCacheFile(new File("test"+"."+"rigid"+".emat"))
                .build()
                .calcEnergyMatrix();

        RCs rcs = sequence.makeRCs(confSpace);

        MultiSequenceSHARKStarBound_refactor.RigidEmatFactory customFact = (SimpleConfSpace customConfSpace) ->{
            // how should we compute energies of molecules?
            EnergyCalculator ecalcMinimized_temp = new EnergyCalculator.Builder(customConfSpace, ffparams)
                    .setParallelism(parallelism)
                    .build();
            // how should we define energies of conformations?
            ConfEnergyCalculator confEcalcMinimized_temp = new ConfEnergyCalculator.Builder(customConfSpace, ecalcMinimized_temp)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(customConfSpace, ecalcMinimized_temp)
                            .build()
                            .calcReferenceEnergies()
                    )
                    .build();

            // BBK* needs rigid energies too
            EnergyCalculator ecalcRigid_temp = new EnergyCalculator.SharedBuilder(ecalcMinimized_temp)
                    .setIsMinimizing(false)
                    .build();
            ConfEnergyCalculator confEcalcRigid_temp = new ConfEnergyCalculator(confEcalcMinimized_temp, ecalcRigid_temp);
            EnergyMatrix rigidEmat_temp = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid_temp)
                    .build()
                    .calcEnergyMatrix();
            return rigidEmat_temp;
        };

        return new MultiSequenceSHARKStarBound_refactor(confSpace, rigidEmat, minEmat, confEcalcMinimized, rcs, parallelism, customFact );

    }

    /**
     * Creates a flexible only confspace with one chain using a few residues from 1CC8
     */
    public SimpleConfSpace make1CC8Flexible(){
        Strand strand1 = new Strand.Builder(metallochaperone).setResidues("A2", "A10").build();
        strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType).setContinuous();

        return new SimpleConfSpace.Builder()
                .addStrands(strand1)
                .build();

    }

    /**
     * Creates a mutable confspace with one chain using 10 residues from 1CC8
     */
    private SimpleConfSpace make1CC8Mutable(){
        Strand strand1 = new Strand.Builder(metallochaperone).setResidues("A2", "A20").build();
        strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType, "ARG").setContinuous();
        strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "ILE").setContinuous();
        strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType, "MET").setContinuous();
        strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A7").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A8").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A9").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A10").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A11").setLibraryRotamers(Strand.WildType).setContinuous();

        return new SimpleConfSpace.Builder()
                .addStrands(strand1)
                .build();

    }

    private SimpleConfSpace make1CC8MutableContinuousSmall() {
        Strand strand1 = new Strand.Builder(metallochaperone).setResidues("A2", "A10").build();
        strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "ILE").setContinuous();
        strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType).setContinuous();

        return new SimpleConfSpace.Builder()
                .addStrands(strand1)
                .build();

    }

    @Test
    public void testMakePfunc(){
        SimpleConfSpace confSpace = make1CC8Flexible();
        Sequence wildType = confSpace.makeWildTypeSequence();
        MultiSequenceSHARKStarBound_refactor bound = makeSHARKStarPfuncForConfSpace(confSpace, wildType, 0.99, null, null);
    }

    @Test
    public void testPrecomputeFlexible(){
        SimpleConfSpace confSpace = make1CC8Mutable();
        Sequence wildType = confSpace.makeWildTypeSequence();
        MultiSequenceSHARKStarBound_refactor bound = makeSHARKStarPfuncForConfSpace(confSpace, wildType, 0.99, null, null);
        bound.init(0.90);
    }

    @Test
    public void testPrecomputeFlexibleTraditional(){
        SimpleConfSpace confSpace = make1CC8Mutable();
        SimpleConfSpace flexCopy = confSpace.makeFlexibleCopy();
        Sequence wildType = flexCopy.makeWildTypeSequence();
        PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(flexCopy, wildType, .68);
        traditionalPfunc.setReportProgress(true);
        traditionalPfunc.compute();

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

        System.out.println(String.format("Gradient Descent pfunc: [%.3f, %.3f]",
                bc.freeEnergy(traditionalPfunc.getValues().calcLowerBound()),
                bc.freeEnergy(traditionalPfunc.getValues().calcUpperBound())));

        traditionalPfunc.printStats();
    }

    @Test
    public void testPrecomputeFlexible_small(){
        SimpleConfSpace confSpace = make1CC8MutableContinuousSmall();
        Sequence wildType = confSpace.makeWildTypeSequence();
        MultiSequenceSHARKStarBound_refactor bound = makeSHARKStarPfuncForConfSpace(confSpace, wildType, 0.68, null, null);
        bound.init(0.68);
    }

    @Test
    public void testComputeForSequence(){
        SimpleConfSpace confSpace = make1CC8Mutable();
        Sequence wildType = confSpace.makeWildTypeSequence();
        MultiSequenceSHARKStarBound_refactor bound = makeSHARKStarPfuncForConfSpace(confSpace, wildType, 0.90, null, null);
        bound.init(0.90);
        PartitionFunction ssbound = bound.getPartitionFunctionForSequence(wildType);
        ssbound.compute();
    }

    @Test
    public void testComputeForSequenceCorrectness(){
        SimpleConfSpace confSpace = make1CC8Mutable();
        Sequence wildType = confSpace.makeWildTypeSequence();

        PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(confSpace, wildType, .1);
        traditionalPfunc.setReportProgress(true);
        traditionalPfunc.compute();

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

        System.out.println(String.format("Gradient Descent pfunc: [%.3f, %.3f]",
                bc.freeEnergy(traditionalPfunc.getValues().calcLowerBound()),
                bc.freeEnergy(traditionalPfunc.getValues().calcUpperBound())));

        traditionalPfunc.printStats();
        //MultiSequenceSHARKStarBound_refactor bound = makeSHARKStarPfuncForConfSpace(confSpace, wildType, 0.99, null, null);
        //bound.init(0.99);
        //PartitionFunction ssbound = bound.getPartitionFunctionForSequence(wildType);
        //ssbound.compute();
    }
}
