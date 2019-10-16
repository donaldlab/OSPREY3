package edu.duke.cs.osprey.sharkstar;

import static edu.duke.cs.osprey.sharkstar.TestSHARKStar.loadFromCFS;
import static edu.duke.cs.osprey.sharkstar.tools.MultiSequenceSHARKStarNodeStatistics.setSigFigs;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.markstar.prototype.SimpleConf;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.jetbrains.annotations.NotNull;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.List;

public class TestSHARKStarBound extends TestBase {

    private static Molecule metallochaperone;
    private static Molecule protein_1a0r;
    private static PartitionFunctionFactory MSSHARKStarPfuncFactory;

    @BeforeClass
    public static void beforeClass() {
        metallochaperone = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");
        //protein_1a0r = PDBIO.readFile("examples//1a0r_prepped.pdb");
    }

    /**
     * Computes a single partition function using <confSpace> for <sequence> to <epsilon>
     */
    private SHARKStarBound makeSHARKStarPfuncForConfSpace(SimpleConfSpace confSpace, @NotNull Sequence sequence, double epsilon, SHARKStarBound preComputedFlex, UpdatingEnergyMatrix preCompCorrections){
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

        // BBK* needs rigid energies too
        EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                .setIsMinimizing(false)
                .build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, ecalcRigid);
        PartitionFunctionFactory pfuncFactory = new PartitionFunctionFactory(confSpace, confEcalcMinimized, "pfunc");
        if (preComputedFlex == null)
            pfuncFactory.setUseSHARKStar(confEcalcRigid);
        else
            pfuncFactory.setUseSHARKStar(confEcalcRigid, preComputedFlex);
        if (preCompCorrections != null)
            pfuncFactory.setPrecomputedCorrections(preCompCorrections);
        // filter the global sequence to this conf space
        // make the partition function
        RCs rcs = sequence.makeRCs(confSpace);

        return (SHARKStarBound) pfuncFactory.makePartitionFunctionFor(rcs, rcs.getNumConformations(), epsilon);

    }

    private PartitionFunction makeMultiSequenceSHARKStarPfuncForConfSpace(SimpleConfSpace confSpace,
                                                          RCs rcs, double epsilon,
                                                          MultiSequenceSHARKStarBound preComputedFlex){
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

            // BBK* needs rigid energies too
            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                    .setIsMinimizing(false)
                    .build();
            ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, ecalcRigid);
            MSSHARKStarPfuncFactory = new PartitionFunctionFactory(confSpace, confEcalcMinimized, "pfunc");
            if (preComputedFlex == null)
                MSSHARKStarPfuncFactory.setUseMSSHARKStar(confEcalcRigid);
            else
                MSSHARKStarPfuncFactory.setUseMSSHARKStar(confEcalcRigid, preComputedFlex);
            // filter the global sequence to this conf space
            // make the partition function
            MSSHARKStarPfuncFactory.setCachePattern("defualt");

            return MSSHARKStarPfuncFactory.makePartitionFunctionFor(rcs, rcs.getNumConformations(), epsilon);
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

    /**
     * Creates a flexible only confspace with one chain using a few residues from 1CC8
     */
    private SimpleConfSpace make1CC8Flexible(){
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

    /**
     * Creates a mutable confspace with one chain using a few residues from 1CC8
     */
    private SimpleConfSpace make1CC8Mutable(){
        Strand strand1 = new Strand.Builder(metallochaperone).setResidues("A2", "A10").build();
        strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
        strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "ILE");
        strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType, "ALA");
        strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType, "MET");
        strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType);
        strand1.flexibility.get("A7").setLibraryRotamers(Strand.WildType);
        strand1.flexibility.get("A8").setLibraryRotamers(Strand.WildType);

        return new SimpleConfSpace.Builder()
                .addStrands(strand1)
                .setShellDistance(9)
                .build();

    }

    /**
     * Creates a mutable confspace with one chain using 10 residues from 1CC8
     */
    private SimpleConfSpace make1CC8MutableContinuous(){
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
                .setShellDistance(9)
                .build();

    }

    private SimpleConfSpace make1CC8MutableContinuousSmall() {
        Strand strand1 = new Strand.Builder(metallochaperone).setResidues("A2", "A10").build();
        strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "ILE").setContinuous();
        strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType).setContinuous();

        return new SimpleConfSpace.Builder()
                .addStrands(strand1)
                .setShellDistance(9)
                .build();

    }

    /**
     * Test that the SHARKStarBound can compute a partition function on a single sequence
     */
    @Test
    public void testComputeSingleSequencePfunc(){
        SimpleConfSpace flexConfSpace = make1CC8Flexible();
        Sequence wildType = flexConfSpace.makeWildTypeSequence();
        PartitionFunction pfunc = makeSHARKStarPfuncForConfSpace(flexConfSpace, wildType, 0.68, null, null);
        pfunc.compute();

        assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
    }

    /**
     * Test that we can make a SHARKStarBound with a precomputed flexible SHARKStarBound
     *
     * Test that the upper and lower pfunc bounds are passed through properly
     */
    @Test
    public void testPrecomputeFlexible(){
        // make full confspace and the flexible copy
        SimpleConfSpace mutableConfSpace = make1CC8Mutable();
        SimpleConfSpace flexCopyConfSpace = mutableConfSpace.makeFlexibleCopy();

        // precompute flexible residues
        Sequence flexSeq = flexCopyConfSpace.makeWildTypeSequence();
        SHARKStarBound preCompFlex = makeSHARKStarPfuncForConfSpace(flexCopyConfSpace, flexSeq, 0.68, null, null);
        preCompFlex.compute();
        BigDecimal precomputedUpper = preCompFlex.getUpperBound();
        BigDecimal precomputedLower = preCompFlex.getLowerBound();

        // make the full confspace partitionFunction
        Sequence fullSeq = mutableConfSpace.makeWildTypeSequence();
        SHARKStarBound fullPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, 0.68, preCompFlex, preCompFlex.genCorrectionMatrix());

        // Test that the precomputed bounds are the same
        assertThat(fullPfunc.getPrecomputedUpperBound(), is(precomputedUpper));
        assertThat(fullPfunc.getPrecomputedLowerBound(), is(precomputedLower));

        // Ensure that the new rootNode bounds are not the same as the precomputed bounds
        assertThat(fullPfunc.getUpperBound(), not(is(fullPfunc.getPrecomputedUpperBound())));
        assertThat(fullPfunc.getLowerBound(), not(is(fullPfunc.getPrecomputedLowerBound())));

    }

    /**
     * Test that we can effectively map the precomputed flexible tree to the full confspace
     */
    @Test
    public void testMappingPrecomputationToFullTree(){
        // make full confspace and the flexible copy
        SimpleConfSpace mutableConfSpace = make1CC8Mutable();
        SimpleConfSpace flexCopyConfSpace = mutableConfSpace.makeFlexibleCopy();

        // precompute flexible residues
        Sequence flexSeq = flexCopyConfSpace.makeWildTypeSequence();
        SHARKStarBound preCompFlex = makeSHARKStarPfuncForConfSpace(flexCopyConfSpace, flexSeq, 0.68, null, null);
        preCompFlex.compute();

        // make the full confspace partitionFunction
        Sequence fullSeq = mutableConfSpace.makeWildTypeSequence();
        SHARKStarBound fullPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, 0.68, preCompFlex, preCompFlex.genCorrectionMatrix());

        int[] expectedArray = {2,4};
        assertThat(fullPfunc.genConfSpaceMapping(), is(expectedArray));
    }

    /**
     * Test to make sure that our epsilon gets reset correctly upon making a precomputed tree
     */
    @Test
    public void testPrecomputedPfuncEpsilon(){
        // make full confspace and the flexible copy
        SimpleConfSpace mutableConfSpace = make1CC8Mutable();
        SimpleConfSpace flexCopyConfSpace = mutableConfSpace.makeFlexibleCopy();

        // precompute flexible residues
        Sequence flexSeq = flexCopyConfSpace.makeWildTypeSequence();
        SHARKStarBound preCompFlex = makeSHARKStarPfuncForConfSpace(flexCopyConfSpace, flexSeq, 0.68, null, null);
        preCompFlex.compute();

        // make the full confspace partitionFunction
        Sequence fullSeq = mutableConfSpace.makeWildTypeSequence();
        SHARKStarBound fullPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, 0.68, preCompFlex, preCompFlex.genCorrectionMatrix());

        assertThat(preCompFlex.epsilonBound, lessThan(0.68));
        assertThat(fullPfunc.epsilonBound, greaterThan(0.68));
    }

    /**
     * Test to make sure that we can compute partition functions after having passed in a flexible precomputed space
     */
    @Test
    public void testPrecomputedPfuncComputation(){
        // make full confspace and the flexible copy
        SimpleConfSpace mutableConfSpace = make1CC8Mutable();
        SimpleConfSpace flexCopyConfSpace = mutableConfSpace.makeFlexibleCopy();

        // precompute flexible residues
        Sequence flexSeq = flexCopyConfSpace.makeWildTypeSequence();
        SHARKStarBound preCompFlex = makeSHARKStarPfuncForConfSpace(flexCopyConfSpace, flexSeq, 0.68, null, null);
        preCompFlex.compute();

        // make the full confspace partitionFunction
        Sequence fullSeq = mutableConfSpace.makeWildTypeSequence();
        SHARKStarBound fullPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, 0.68, preCompFlex, preCompFlex.genCorrectionMatrix());

        // update and compute
        fullPfunc.compute();
    }

    /**
     * Test to make sure that the 2-step computed partition functions return the same bounds as the one-step pfunc
     */
    @Test
    public void testPreComputedPfuncCorrectness(){
        double epsilon = 0.68;
        // make full confspace and the flexible copy
        SimpleConfSpace mutableConfSpace = make1CC8Mutable();
        SimpleConfSpace flexCopyConfSpace = mutableConfSpace.makeFlexibleCopy();

        // precompute flexible residues
        Sequence flexSeq = flexCopyConfSpace.makeWildTypeSequence();
        SHARKStarBound preCompFlex = makeSHARKStarPfuncForConfSpace(flexCopyConfSpace, flexSeq, epsilon, null, null);
        preCompFlex.compute();

        // make the full confspace partitionFunction, and compute it much, much more accurately.
        Sequence fullSeq = mutableConfSpace.makeWildTypeSequence();
        SHARKStarBound fullPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, epsilon, preCompFlex, preCompFlex.genCorrectionMatrix());

        // update and compute
        //fullPfunc.updatePrecomputedConfTree();
        fullPfunc.compute();

        PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(mutableConfSpace, fullSeq, epsilon/10000);
        traditionalPfunc.setReportProgress(true);
        traditionalPfunc.compute();

        System.out.println("=============================================================================================");
        System.out.println("====================== Running SHARK* without precomputed tree ==============================");
        System.out.println("=============================================================================================");
        // compute partition function the regular way
        SHARKStarBound regPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, epsilon, null, null);
        regPfunc.compute();

        System.out.println("Gradient Descent pfunc: "+formatBounds(traditionalPfunc.getValues().calcLowerBound(),
                traditionalPfunc.getValues().calcUpperBound()));
        System.out.println("RegPfunc: " + formatBounds(regPfunc.getValues().calcLowerBound(), regPfunc.getValues().calcUpperBound()));
        System.out.println("precompPfunc: " + formatBounds(fullPfunc.getValues().calcLowerBound(), fullPfunc.getValues().calcUpperBound()));
        assertThat(fullPfunc.getStatus(), is(PartitionFunction.Status.Estimated));
        assertThat(regPfunc.getStatus(), is(PartitionFunction.Status.Estimated));
        assertThat(regPfunc.getValues().calcUpperBound(), greaterThan(traditionalPfunc.getValues().calcUpperBound()));
        assertThat(regPfunc.getValues().calcLowerBound(), lessThan(traditionalPfunc.getValues().calcLowerBound()));
        assertThat(fullPfunc.getValues().calcUpperBound(), greaterThan(traditionalPfunc.getValues().calcUpperBound()));
        assertThat(fullPfunc.getValues().calcLowerBound(), lessThan(traditionalPfunc.getValues().calcLowerBound()));

    }

    @Test
    public void testMultiSequencePFuncFactory () {
        double epsilon = 0.01;
        // make full confspace and the flexible copy
        SimpleConfSpace mutableConfSpace = make1CC8Mutable();

        Sequence fullSeq = mutableConfSpace.makeUnassignedSequence();
        MultiSequenceSHARKStarBound fullPfunc =
                (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(mutableConfSpace,
                        fullSeq.makeRCs(mutableConfSpace), epsilon, null);

        PartitionFunction wtBound =
                fullPfunc.getPartitionFunctionForSequence(mutableConfSpace.makeWildTypeSequence());
        wtBound.compute();

        System.out.println("========================== Now computing mutant sequence ========================");
        Sequence mutantSequence = mutableConfSpace.makeWildTypeSequence() .set("A5","MET") .set("A3","ILE");
        RCs mutantRcs = mutantSequence.makeRCs(mutableConfSpace);
        PartitionFunction muttBound =
                MSSHARKStarPfuncFactory.makePartitionFunctionFor(mutantRcs, mutantRcs.getNumConformations(), epsilon, mutantSequence);
        muttBound.compute();

        PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(mutableConfSpace, mutantSequence, epsilon);
        traditionalPfunc.setReportProgress(true);
        traditionalPfunc.compute();
        System.out.println("Gradient Descent pfunc: "+formatBounds(traditionalPfunc.getValues().calcLowerBound(),
                traditionalPfunc.getValues().calcUpperBound()));
        System.out.println("precompPfunc: " + formatBounds(muttBound.getValues().calcLowerBound(), muttBound.getValues().calcUpperBound()));

        assertThat(wtBound.getStatus(), is(PartitionFunction.Status.Estimated));
        assertThat(muttBound.getStatus(), is(PartitionFunction.Status.Estimated));
    }

    @Test
    public void testMultiSequence() {
        double epsilon = 0.01;
        // make full confspace and the flexible copy
        SimpleConfSpace mutableConfSpace = make1CC8Mutable();
/*
        SimpleConfSpace flexCopyConfSpace = mutableConfSpace.makeFlexibleCopy();

        // precompute flexible residues
        Sequence flexSeq = flexCopyConfSpace.makeWildTypeSequence();
        MultiSequenceSHARKStarBound preCompFlex = (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(
                flexCopyConfSpace, flexSeq.makeRCs(flexCopyConfSpace), epsilon, null);
        PartitionFunction ssbound = preCompFlex.getPartitionFunctionForSequence(flexSeq);
        ssbound.compute();

        System.out.println("Precomputed flex done.");

 */
        // make the full confspace partitionFunction, and compute it much, much more accurately.
        Sequence fullSeq = mutableConfSpace.makeUnassignedSequence();
        MultiSequenceSHARKStarBound fullPfunc =
                (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(mutableConfSpace,
                        fullSeq.makeRCs(mutableConfSpace), epsilon, null);

        PartitionFunction wtBound =
                fullPfunc.getPartitionFunctionForSequence(mutableConfSpace.makeWildTypeSequence());
        wtBound.compute();

        System.out.println("========================== Now computing mutant sequence ========================");
        Sequence mutantSequence = mutableConfSpace.makeWildTypeSequence() .set("A5","MET") .set("A3","ILE");
        PartitionFunction muttBound =
                fullPfunc.getPartitionFunctionForSequence(mutantSequence);
        muttBound.compute();

        PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(mutableConfSpace, mutantSequence, epsilon);
        traditionalPfunc.setReportProgress(true);
        traditionalPfunc.compute();
        System.out.println("Gradient Descent pfunc: "+formatBounds(traditionalPfunc.getValues().calcLowerBound(),
                traditionalPfunc.getValues().calcUpperBound()));
        System.out.println("precompPfunc: " + formatBounds(muttBound.getValues().calcLowerBound(), muttBound.getValues().calcUpperBound()));

        assertThat(wtBound.getStatus(), is(PartitionFunction.Status.Estimated));
        assertThat(muttBound.getStatus(), is(PartitionFunction.Status.Estimated));
    }

    @Test
    public void testContinuous_3ma2() {

        double epsilon = 0.68;
        try {
            SimpleConfSpace mutableConfSpace = loadFromCFS("test-resources/3ma2_A_6res_3.157E+06.cfs").complex;
            Sequence fullSeq = mutableConfSpace.makeUnassignedSequence();
            MultiSequenceSHARKStarBound fullPfunc =
                    (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(mutableConfSpace,
                            fullSeq.makeRCs(mutableConfSpace), epsilon, null);

            PartitionFunction wtBound =
                    fullPfunc.getPartitionFunctionForSequence(mutableConfSpace.makeWildTypeSequence());
            wtBound.compute();

            System.out.println("========================== Now computing mutant sequence ========================");
            Sequence mutantSequence = mutableConfSpace.makeWildTypeSequence() .set("A5","MET") .set("A3","ILE");
            PartitionFunction muttBound =
                    fullPfunc.getPartitionFunctionForSequence(mutantSequence);
            muttBound.compute();

            PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(mutableConfSpace, mutantSequence, epsilon);
            traditionalPfunc.setReportProgress(true);
            traditionalPfunc.compute();
            System.out.println("Gradient Descent pfunc: "+formatBounds(traditionalPfunc.getValues().calcLowerBound(),
                    traditionalPfunc.getValues().calcUpperBound()));
            System.out.println("precompPfunc: " + formatBounds(muttBound.getValues().calcLowerBound(), muttBound.getValues().calcUpperBound()));
        }catch(IOException e) {

        }
    }

    @Test
    public void test1CC8Small() {
        double epsilon = 0.68;
        SimpleConfSpace mutableConfSpace = make1CC8MutableContinuousSmall();
        Sequence fullSeq = mutableConfSpace.makeUnassignedSequence();
        MultiSequenceSHARKStarBound fullPfunc =
                (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(mutableConfSpace,
                        fullSeq.makeRCs(mutableConfSpace), epsilon, null);

        PartitionFunction wtBound =
                fullPfunc.getPartitionFunctionForSequence(mutableConfSpace.makeWildTypeSequence());
        wtBound.compute();

        System.out.println("========================== Now computing mutant sequence ========================");
        Sequence mutantSequence = mutableConfSpace.makeWildTypeSequence() .set("A3","ILE");
        PartitionFunction muttBound =
                fullPfunc.getPartitionFunctionForSequence(mutantSequence);
        muttBound.compute();

        PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(mutableConfSpace, mutantSequence, epsilon);
        traditionalPfunc.setReportProgress(true);
        traditionalPfunc.compute();
        System.out.println("Gradient Descent pfunc: "+formatBounds(traditionalPfunc.getValues().calcLowerBound(),
                traditionalPfunc.getValues().calcUpperBound()));
        System.out.println("precompPfunc: " + formatBounds(muttBound.getValues().calcLowerBound(), muttBound.getValues().calcUpperBound()));
    }

    @Test
    public void testMultiSequenceContinuous() {

        double epsilon = 0.68;
        SimpleConfSpace mutableConfSpace = make1CC8MutableContinuous();
        // make full confspace and the flexible copy
        SimpleConfSpace flexCopyConfSpace = mutableConfSpace.makeFlexibleCopy();

        // precompute flexible residues
        Sequence flexSeq = flexCopyConfSpace.makeWildTypeSequence();
        MultiSequenceSHARKStarBound preCompFlex = (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(
                flexCopyConfSpace, flexSeq.makeRCs(flexCopyConfSpace), epsilon, null);
        PartitionFunction ssbound = preCompFlex.getPartitionFunctionForSequence(flexSeq);
        ssbound.compute();

        System.out.println("Precomputed flex done.");

        // make the full confspace partitionFunction, and compute it much, much more accurately.
        Sequence fullSeq = mutableConfSpace.makeUnassignedSequence();
        MultiSequenceSHARKStarBound fullPfunc =
                (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(mutableConfSpace,
                        fullSeq.makeRCs(mutableConfSpace), epsilon, preCompFlex);

        System.out.println("========================== Now computing wildtype sequence ========================");
        PartitionFunction wtBound =
                fullPfunc.getPartitionFunctionForSequence(mutableConfSpace.makeWildTypeSequence());
        wtBound.compute();

        System.out.println("========================== Now computing mutant sequence ========================");
        Sequence mutantSequence = mutableConfSpace.makeWildTypeSequence() .set("213","LYS");
        PartitionFunction muttBound =
                fullPfunc.getPartitionFunctionForSequence(mutantSequence);
        muttBound.compute();

        PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(mutableConfSpace, mutantSequence, epsilon);
        traditionalPfunc.setReportProgress(true);
        traditionalPfunc.compute();
        System.out.println("Gradient Descent pfunc: "+formatBounds(traditionalPfunc.getValues().calcLowerBound(),
                traditionalPfunc.getValues().calcUpperBound()));
        System.out.println("precompPfunc: " + formatBounds(muttBound.getValues().calcLowerBound(), muttBound.getValues().calcUpperBound()));

        assertThat(wtBound.getStatus(), is(PartitionFunction.Status.Estimated));
        assertThat(muttBound.getStatus(), is(PartitionFunction.Status.Estimated));

    }

    @Test
    public void debug3ma2Bigger() {
        double epsilon = 0.68;
        try {
            SimpleConfSpace mutableConfSpace = loadFromCFS("test-resources/3ma2_D_10res_1.140E+10.cfs").complex;
            Sequence fullSeq = mutableConfSpace.makeUnassignedSequence();
            MultiSequenceSHARKStarBound fullPfunc =
                    (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(mutableConfSpace,
                            fullSeq.makeRCs(mutableConfSpace), epsilon, null);

            boolean runWT = true;
            if(runWT) {
                PartitionFunction wtBound =
                        fullPfunc.getPartitionFunctionForSequence(mutableConfSpace.makeWildTypeSequence());
                wtBound.compute();
            }

            System.out.println("========================== Now computing mutant sequence ========================");
            Sequence mutantSequence = mutableConfSpace.makeWildTypeSequence() .set("D198","GLY");
            PartitionFunction muttBound =
                    fullPfunc.getPartitionFunctionForSequence(mutantSequence);
            muttBound.compute();

            PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(mutableConfSpace, mutantSequence, epsilon);
            traditionalPfunc.setReportProgress(true);
            traditionalPfunc.compute();
            System.out.println("Gradient Descent pfunc: "+formatBounds(traditionalPfunc.getValues().calcLowerBound(),
                    traditionalPfunc.getValues().calcUpperBound()));
            System.out.println("precompPfunc: " + formatBounds(muttBound.getValues().calcLowerBound(), muttBound.getValues().calcUpperBound()));
            }
        catch(IOException e) {
                e.printStackTrace();

        }
    }

    @Test
    public void debug3u7y() {
        double epsilon = 0.68;
        try {
            SimpleConfSpace mutableConfSpace = loadFromCFS("test-resources/3u7y_L_4res_9.135E+04.cfs").complex;
            Sequence fullSeq = mutableConfSpace.makeUnassignedSequence();
            MultiSequenceSHARKStarBound fullPfunc =
                    (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(mutableConfSpace,
                            fullSeq.makeRCs(mutableConfSpace), epsilon, null);

            boolean runWT = true;
            if(runWT) {
                PartitionFunction wtBound =
                        fullPfunc.getPartitionFunctionForSequence(mutableConfSpace.makeWildTypeSequence());
                wtBound.compute();
            }

            System.out.println("========================== Now computing mutant sequence ========================");
            Sequence mutantSequence = mutableConfSpace.makeWildTypeSequence() .set("L117","LEU");
            PartitionFunction muttBound =
                    fullPfunc.getPartitionFunctionForSequence(mutantSequence);
            muttBound.compute();

            PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(mutableConfSpace, mutantSequence, epsilon);
            traditionalPfunc.setReportProgress(true);
            traditionalPfunc.compute();
            System.out.println("Gradient Descent pfunc: "+formatBounds(traditionalPfunc.getValues().calcLowerBound(),
                    traditionalPfunc.getValues().calcUpperBound()));
            System.out.println("precompPfunc: " + formatBounds(muttBound.getValues().calcLowerBound(), muttBound.getValues().calcUpperBound()));
        }
        catch(IOException e) {
            e.printStackTrace();

        }
    }

    private SimpleConfSpace make1CC8MutableContinuous_debug() {
    Strand strand1 = new Strand.Builder(metallochaperone).setResidues("A2", "A10").build();
        strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType, "ALA", "ARG").setContinuous();
        strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "GLU", "ILE").setContinuous();
        strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType, "LYS", "MET").setContinuous();
        strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType).setContinuous();

        return new SimpleConfSpace.Builder()
                .addStrands(strand1)
                .setShellDistance(9)
                .build();

    }

    private SimpleConfSpace make1gua_debug() {
        return TestKStar.make1GUA11().complex;
    }

    private SimpleConfSpace make1a0r_debug() {
        Strand strand1 = new Strand.Builder(protein_1a0r).setResidues("B2", "B340").build();
        strand1.flexibility.get("B282").setLibraryRotamers(Strand.WildType, "ALA", "ARG").addWildTypeRotamers().setContinuous();
        strand1.flexibility.get("B235").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        return new SimpleConfSpace.Builder()
                .addStrands(strand1)
                .setShellDistance(9)
                .build();
    }

    public void testMultiSequenceCorrectness() {
        double epsilon = 0.68;
        // make full confspace and the flexible copy
        SimpleConfSpace mutableConfSpace = make1CC8Mutable();
        SimpleConfSpace flexCopyConfSpace = mutableConfSpace.makeFlexibleCopy();

        // precompute flexible residues
        Sequence flexSeq = flexCopyConfSpace.makeWildTypeSequence();
        SHARKStarBound preCompFlex = makeSHARKStarPfuncForConfSpace(flexCopyConfSpace, flexSeq, epsilon, null, null);
        preCompFlex.compute();

        // make the full confspace partitionFunction, and compute it much, much more accurately.
        Sequence fullSeq = mutableConfSpace.makeWildTypeSequence();
        SHARKStarBound fullPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, epsilon, preCompFlex, preCompFlex.genCorrectionMatrix());

        // update and compute
        //fullPfunc.updatePrecomputedConfTree();
        fullPfunc.compute();

        PartitionFunction traditionalPfunc = makeGradientDescentPfuncForConfSpace(mutableConfSpace, fullSeq, epsilon/10000);
        traditionalPfunc.setReportProgress(true);
        traditionalPfunc.compute();

        System.out.println("=============================================================================================");
        System.out.println("====================== Running SHARK* without precomputed tree ==============================");
        System.out.println("=============================================================================================");
        // compute partition function the regular way
        SHARKStarBound regPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, epsilon, null, null);
        regPfunc.compute();

        System.out.println("Gradient Descent pfunc: "+formatBounds(traditionalPfunc.getValues().calcLowerBound(),
                traditionalPfunc.getValues().calcUpperBound()));
        System.out.println("RegPfunc: " + formatBounds(regPfunc.getValues().calcLowerBound(), regPfunc.getValues().calcUpperBound()));
        System.out.println("precompPfunc: " + formatBounds(fullPfunc.getValues().calcLowerBound(), fullPfunc.getValues().calcUpperBound()));
        assertThat(fullPfunc.getStatus(), is(PartitionFunction.Status.Estimated));
        assertThat(regPfunc.getStatus(), is(PartitionFunction.Status.Estimated));
        assertThat(regPfunc.getValues().calcUpperBound(), greaterThan(traditionalPfunc.getValues().calcUpperBound()));
        assertThat(regPfunc.getValues().calcLowerBound(), lessThan(traditionalPfunc.getValues().calcLowerBound()));
        assertThat(fullPfunc.getValues().calcUpperBound(), greaterThan(traditionalPfunc.getValues().calcUpperBound()));
        assertThat(fullPfunc.getValues().calcLowerBound(), lessThan(traditionalPfunc.getValues().calcLowerBound()));

    }

    private String formatBounds(BigDecimal lower, BigDecimal upper) {
        return "[" + setSigFigs(lower) + "," + setSigFigs(upper) + "]";
    }

    @Test
    public void testCapturingPartialMinimizations(){
        double epsilon = 0.68;
        // make full confspace and the flexible copy
        SimpleConfSpace mutableConfSpace = make1CC8MutableContinuous();
        SimpleConfSpace flexCopyConfSpace = mutableConfSpace.makeFlexibleCopy();

        // precompute flexible residues
        Sequence flexSeq = flexCopyConfSpace.makeWildTypeSequence();
        SHARKStarBound preCompFlex = makeSHARKStarPfuncForConfSpace(flexCopyConfSpace, flexSeq, epsilon, null, null);
        preCompFlex.compute();

        // get the list of corrections
        List<TupE> flexibleCorrections = preCompFlex.correctionMatrix.getAllCorrections();

        UpdatingEnergyMatrix precomputedCorrections = preCompFlex.genCorrectionMatrix();
        // Pre-fullpfunc minlist
        System.out.println(preCompFlex.minList.toString());

        // make the full confspace partitionFunction, and compute it much, much more accurately.
        Sequence fullSeq = mutableConfSpace.makeWildTypeSequence();
        SHARKStarBound fullPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, epsilon, preCompFlex, precomputedCorrections);

        // update and compute
        fullPfunc.compute();

        // Check that the number of corrections is right
        assertThat((precomputedCorrections.getTrieSize()), is(fullPfunc.TestNumCorrections));

        // Check that all of the corrections are properly stored in the new list
        RCTuple permTup = null;
        for (TupE tupe : flexibleCorrections){
            // assert that all of them exist
            permTup = tupe.tup.permutedCopy(fullPfunc.getConfSpacePermutation());
            assertThat(fullPfunc.correctionMatrix.hasHigherOrderTermFor(permTup), is(true));
            System.out.println(permTup);
            // assert that the energies are correct
            System.out.println(fullPfunc.correctionMatrix.getInternalEnergy(permTup));
        }
    }


}

