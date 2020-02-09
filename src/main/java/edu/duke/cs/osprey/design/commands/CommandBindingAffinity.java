package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.models.AffinityDesign;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

@Parameters(commandDescription = CommandBindingAffinity.CommandDescription)
public class CommandBindingAffinity extends RunnableCommand {

    public static final String CommandName = "affinity";
    public static final String CommandDescription = "Compute an epsilon approximation to binding affinity (K*).";

    @Override
    public int run(JCommander commander, String[] args) {
        var opt = processHelpAndNoArgs(commander, args);
        if (opt.isPresent()) {
            return opt.get();
        }

        return parseAndValidate(delegate.design)
                .map(this::runAffinityDesign)
                .orElse(Main.Failure);
    }

    private int runAffinityDesign(AffinityDesign design) {
        var forcefieldParams = new ForcefieldParams();
        var confSpace1 = delegate.createConfSpace(design.protein, forcefieldParams);
        var confSpace2 = delegate.createConfSpace(design.ligand, forcefieldParams);
        var strands = List.of(confSpace1.strands, confSpace2.strands)
                .stream()
                .flatMap(List::stream)
                .collect(Collectors.toList());

        var complexConfSpace = new SimpleConfSpace.Builder()
                .addStrands(strands)
                .build();

        // Exit early if just trying to validate input
        if (delegate.verifyInput) {
            System.out.println("Design file validated.");
            return Main.Success;
        }

        /* Decides whether to use CPU(s) and/or GPU(s) (purely implementation specific) */
        var parallelism = delegate.getParallelism();

        /* Used to calculate energies of a molecule, also used to minimize the molecule */
        var energyCalculator = new EnergyCalculator.Builder(complexConfSpace, forcefieldParams)
                .setParallelism(parallelism)
                .build();

        var epsilon = delegate.epsilon > 0 ? delegate.epsilon : design.epsilon;
        var kstar = new KStar(confSpace1, confSpace2, complexConfSpace, makeKStarSettings(epsilon));

        for (var info : kstar.confSpaceInfos()) {
            var referenceEnergies = new SimpleReferenceEnergies.Builder(info.confSpace, energyCalculator).build();

            info.confEcalc = new ConfEnergyCalculator.Builder(info.confSpace, energyCalculator)
                    .setReferenceEnergies(referenceEnergies)
                    .build();

            var energyMatrix = new SimplerEnergyMatrixCalculator.Builder(info.confSpace, energyCalculator)
                    .build()
                    .calcEnergyMatrix();

            info.confSearchFactory = rcs -> new ConfAStarTree.Builder(energyMatrix, rcs)
                    .setShowProgress(false)
                    .build();
        }

        printResults(kstar.run());
        return Main.Success;
    }

    private void printResults(List<KStar.ScoredSequence> results) {
        for (var result : results) {
            System.out.println("result: ");
            System.out.println(String.format("\tsequence: %s", result.sequence));
            System.out.println(String.format("\tscore: %s", result.score));
        }
    }

    private KStar.Settings makeKStarSettings(double epsilon) {
        return new KStar.Settings.Builder()
                .setEpsilon(epsilon)
                .addScoreConsoleWriter()
                .build();
    }

    static Optional<AffinityDesign> parseAndValidate(File designSpec) {

        AffinityDesign design;
        try {
            design = AffinityDesign.parse(designSpec);
        } catch (IOException e) {
            e.printStackTrace();
            return Optional.empty();
        }

        var specErrors = design.validate();
        if (!specErrors.isEmpty()) {
            System.err.println("Invalid design specification. The following validations failed:");
            specErrors.stream().map(s -> String.format("- %s", s)).forEach(System.err::println);
            return Optional.empty();
        }

        return Optional.of(design);
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
