package edu.duke.cs.osprey.design.commands;


import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.models.AminoAcid;
import edu.duke.cs.osprey.design.models.StabilityDesign;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.IOException;
import java.util.stream.Collectors;

@Parameters(commandDescription = CommandPartitionFunction.CommandDescription)
public class CommandPartitionFunction extends RunnableCommand {

    public static final String CommandName = "stability";
    static final String CommandDescription = "Estimate the partition function value(s) of different conformations";

    @Override
    public int run(JCommander commander, String[] args) {
        var retVal = processHelpAndNoArgs(commander, args);

        if (retVal.isPresent()) {
            return retVal.get();
        }

        StabilityDesign design;

        try {
            design = StabilityDesign.parse(delegate.design);
        } catch (IOException e) {
            e.printStackTrace();
            return Main.Failure;
        }

        /* Reads a PDB file into a Molecule. */
        var molecule = PDBIO.read(design.molecule);

        /* This reads parm96a.dat, which contains  the energy parameters of DNA, RNA, and protein residues */
        var ffParams = new ForcefieldParams();

        /* Reads the templates, rotamers, and entropy info for a given forcefield */
        /*
            "/config/parm96a.dat",
            "/config/all_amino94.in",
            "/config/all_aminont94.in",
            "/config/all_aminoct94.in",
            "/config/all_nuc94_and_gr.in",
         */
        var templateLibrary = new ResidueTemplateLibrary.Builder(ffParams.forcefld)
                .build();

        /* Strands combine a Molecule with design flexibility and templates */
        var protein = new Strand.Builder(molecule)
                .setTemplateLibrary(templateLibrary)
                .build();

        /* This block makes certain residues in the protein flexible or mutable */
        design.residueModifiers.forEach(mod -> {
            var residue = protein.flexibility.get(Integer.toString(mod.identity.residueNumber));
            residue.addWildTypeRotamers = true; // mod.flexibility.includeStructureRotamer;
            var toMutations = mod.mutability.stream().map(AminoAcid::toValue).collect(Collectors.toUnmodifiableList());

            // null is okay, as is a List, but not an empty List :(
            if (!toMutations.isEmpty()) {
                residue.setLibraryRotamers(toMutations);
            } else {
                residue.setLibraryRotamers();
            }
        });

        /* Maintains flexibility information with the molecule, and can use that to make new molecules */
        var confSpace = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .setShellDistance(0)
                .build();

        /* Decides whether to use CPU(s) and/or GPU(s) (only about runtime specifications) */
        var parallelism = new Parallelism(Runtime.getRuntime().availableProcessors(), 0, 0);

        /* Used to calculate energies of a molecule, also used to minimize the molecule */
        var energyCalculator = new EnergyCalculator.Builder(confSpace, ffParams)
                .setParallelism(parallelism)
                .build();

        /*
         * Calculate energy for molecules created from conformation spaces.
         *
         * Provides support for applying conformation energy modifications,
         * such as reference energies, residue entropies, and energy partitions.
         */
        var confEnergyCalculator = new ConfEnergyCalculator.Builder(confSpace, energyCalculator)
                .build();

        /* Contains the confSpace and a pruning matrix */
        var rcs = new RCs(confSpace);

        var partitionFnBuilder = new PartitionFunctionFactory(confSpace, confEnergyCalculator, "default");
        partitionFnBuilder.setUseGradientDescent();

        var partFn = partitionFnBuilder.makePartitionFunctionFor(rcs, null, 0.3);
        partFn.compute();
        var evaluated = partFn.getNumConfsEvaluated();
        return Main.Success;
    }

    @Override
    public String getCommandName() {
        return CommandName;
    }

    @Override
    public String getCommandDescription() {
        return CommandDescription;
    }
}
