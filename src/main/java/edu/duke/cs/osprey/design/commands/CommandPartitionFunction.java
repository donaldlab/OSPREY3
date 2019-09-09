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

        var molecule = PDBIO.read(design.molecule);
        var ffParams = new ForcefieldParams();
        var templateLibrary = new ResidueTemplateLibrary.Builder(ffParams.forcefld)
                .build();
        var protein = new Strand.Builder(molecule)
                .setTemplateLibrary(templateLibrary)
                .build();

        design.residueModifiers.forEach(mod -> {
            var residue = protein.flexibility.get(Integer.toString(mod.identity.residueNumber));
            residue.addWildTypeRotamers = true; // mod.flexibility.includeStructureRotamer;
            var toMutations = mod.mutability.stream().map(AminoAcid::toValue).collect(Collectors.toUnmodifiableList());

            if (!toMutations.isEmpty()) {
                residue.setLibraryRotamers(toMutations);
            } else {
                residue.setLibraryRotamers();
            }
        });

        var confSpace = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();

        var parallelism = new Parallelism(Runtime.getRuntime().availableProcessors(), 0, 0);

        var energyCalculator = new EnergyCalculator.Builder(confSpace, ffParams)
                .setParallelism(parallelism)
                .build();

        var confEnergyCalculator = new ConfEnergyCalculator.Builder(confSpace, energyCalculator)
                .build();

        var partitionFnBuilder = new PartitionFunctionFactory(confSpace, confEnergyCalculator, "default");
        partitionFnBuilder.setUseGradientDescent();

        var rcs = new RCs(confSpace);
        var partFn = partitionFnBuilder.makePartitionFunctionFor(rcs, null, design.epsilon);
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
