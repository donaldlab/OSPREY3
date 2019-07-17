package edu.duke.cs.osprey.design;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.MissingCommandException;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import edu.duke.cs.osprey.design.commands.CommandBindingAffinity;
import edu.duke.cs.osprey.design.commands.CommandPartitionFunction;

import java.util.Map;

import static java.lang.System.exit;

public class Main {

    @Parameter(names = {"--help", "-h"}, help = true)
    private boolean help;

    public final static int Success = 0;
    public final static int Failure = 1;
    public final static String ProgramName = "osprey";

    private void printErrorMessage(String msg) {
        System.out.println();
        System.out.println(msg);
        System.out.println();
        System.out.println("Use --help or <command> --help for more info.");
    }

    public int run(String[] args) {

        var commandMap = Map.of(
                CommandPartitionFunction.CommandName, new CommandPartitionFunction(),
                CommandBindingAffinity.CommandName, new CommandBindingAffinity()
        );

        var builder = JCommander.newBuilder()
                .programName(ProgramName)
                .addObject(this);
        commandMap.forEach((name, command) -> builder.addCommand(name, command));
        var commander = builder.build();

        try {
            commander.parse(args);
        } catch (MissingCommandException ex) {
            return Failure;
        } catch (ParameterException ex) {
            printErrorMessage(String.format("Error: %s", ex.getMessage()));
            return Failure;
        }

        if (args.length == 0 || this.help) {
            commander.usage();
            return Success;
        }

        var parsedCommand = commander.getParsedCommand();
        return commandMap.get(parsedCommand).run(commander.getCommands().get(parsedCommand), args);
    }

    public static void main(String[] args) {
        exit(new Main().run(args));

        /*
        var designFile = new File(args[0]);
        var stabilityDesign = StabilityDesign.parse(designFile);
        var molecule = PDBIO.read(stabilityDesign.molecule);

        var ffParams = new ForcefieldParams();
        var templateLibrary = new ResidueTemplateLibrary.Builder(ffParams.forcefld)
                .build();

        var protein = new Strand.Builder(molecule)
                .setTemplateLibrary(templateLibrary)
                .build();

        stabilityDesign.residueModifiers.forEach(mod -> {
            var residue = protein.flexibility.get(mod.identity.positionIdentifier());
            residue.addWildTypeRotamers = mod.flexibility.includeStructureRotamer;
            var toMutations = mod.mutable.stream().map(AminoAcid::toValue).collect(Collectors.toUnmodifiableList());
            if (!toMutations.isEmpty()) {
                residue.setLibraryRotamers(toMutations);
            }
            residue.setContinuous();
        });

        var confSpace = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .setShellDistance(4)
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
        var partFn = partitionFnBuilder.makePartitionFunctionFor(rcs, null, stabilityDesign.epsilon);
        partFn.compute();
        var evaluated = partFn.getNumConfsEvaluated();
         */
    }
}

