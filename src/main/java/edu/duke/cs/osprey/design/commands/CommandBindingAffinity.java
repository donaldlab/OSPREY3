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
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;

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
        var minimizingECalc = new EnergyCalculator.Builder(complexConfSpace, forcefieldParams)
                .setParallelism(parallelism)
                .build();

        var rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingECalc)
                .setIsMinimizing(false)
                .build();

        var epsilon = delegate.epsilon > 0 ? delegate.epsilon : design.epsilon;
        var bbkstar = new BBKStar(confSpace1, confSpace2, complexConfSpace, makeKStarSettings(epsilon), makeBBKStarSettings());

        for (var info : bbkstar.confSpaceInfos()) {

            var referenceEnergies = new SimpleReferenceEnergies.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc).build();
            var confECalcMinimized = new ConfEnergyCalculator.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc)
                    .setReferenceEnergies(referenceEnergies)
                    .build();

            var minimizedEnergyMatrix = new SimplerEnergyMatrixCalculator.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc)
                    .build()
                    .calcEnergyMatrix();

            var rigidEnergyMatrix = new SimplerEnergyMatrixCalculator.Builder(((SimpleConfSpace) info.confSpace), rigidEcalc)
                    .build()
                    .calcEnergyMatrix();

            info.confEcalcMinimized = confECalcMinimized;
            info.confSearchFactoryMinimized = rcs -> new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build();
            info.confSearchFactoryRigid = rcs -> new ConfAStarTree.Builder(rigidEnergyMatrix, rcs).setTraditional().build();

            info.pfuncFactory = (rcs) -> new GradientDescentPfunc(
                    info.confEcalcMinimized,
                    new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                    new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                    rcs.getNumConformations()
            );
        }

        printResults(bbkstar.run(minimizingECalc.tasks));
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

    private BBKStar.Settings makeBBKStarSettings() {
        return new BBKStar.Settings.Builder().build();
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
