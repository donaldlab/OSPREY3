package edu.duke.cs.osprey.sharkstar;

import static org.hamcrest.Matchers.*;
import static org.hamcrest.collection.IsIterableContainingInAnyOrder.containsInAnyOrder;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.jetbrains.annotations.NotNull;
import org.junit.BeforeClass;
import org.junit.Test;

public class TestSHARKStarBound extends TestBase {

    private static Molecule metallochaperone;

    @BeforeClass
    public static void beforeClass() {
        metallochaperone = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");
    }

    /**
     * Computes a single partition function using <confSpace> for <sequence> to <epsilon>
     */
    private SHARKStarBound makeSHARKStarPfuncForConfSpace(SimpleConfSpace confSpace, @NotNull Sequence sequence, double epsilon, SHARKStarBound preComputedFlex){
        // Set up partition function requirements
        Parallelism parallelism = Parallelism.makeCpu(4);
        ForcefieldParams ffparams = new ForcefieldParams();

        // how should we compute energies of molecules?
        try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build()) {
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
            // filter the global sequence to this conf space
            // make the partition function
            RCs rcs = sequence.makeRCs(confSpace);

            return (SHARKStarBound) pfuncFactory.makePartitionFunctionFor(rcs, rcs.getNumConformations(), epsilon);
        }
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
        strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType, "ARG");
        strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "ILE");
        strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
        strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType, "MET");
        strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType);

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
        PartitionFunction pfunc = makeSHARKStarPfuncForConfSpace(flexConfSpace, wildType, 0.68, null);
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
        SHARKStarBound preCompFlex = makeSHARKStarPfuncForConfSpace(flexCopyConfSpace, flexSeq, 0.68, null);
        preCompFlex.compute();

        // make the full confspace partitionFunction
        Sequence fullSeq = mutableConfSpace.makeWildTypeSequence();
        SHARKStarBound fullPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, 0.68, preCompFlex);

        // Test that the precomputed bounds are the same
        assertThat(fullPfunc.getPrecomputedUpperBound(), is(preCompFlex.getUpperBound()));
        assertThat(fullPfunc.getPrecomputedLowerBound(), is(preCompFlex.getLowerBound()));

        // Ensure that the new rootNode bounds are not the same
        assertThat(fullPfunc.getUpperBound(), not(is(preCompFlex.getUpperBound())));
        assertThat(fullPfunc.getLowerBound(), not(is(preCompFlex.getLowerBound())));

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
        SHARKStarBound preCompFlex = makeSHARKStarPfuncForConfSpace(flexCopyConfSpace, flexSeq, 0.68, null);
        preCompFlex.compute();

        // make the full confspace partitionFunction
        Sequence fullSeq = mutableConfSpace.makeWildTypeSequence();
        SHARKStarBound fullPfunc = makeSHARKStarPfuncForConfSpace(mutableConfSpace, fullSeq, 0.68, preCompFlex);

        //clear screen sort of
        System.out.println("\n\n\n\n\n\n\n####################");
        // map the precomputed tree onto the full confspace tree
        fullPfunc.updateConfTree();

    }
}

