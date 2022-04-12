package edu.duke.cs.osprey.design.commands;


import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.analysis.CommandAnalysis;
import edu.duke.cs.osprey.design.analysis.EnergyAnalysisConfListener;
import edu.duke.cs.osprey.design.analysis.ThermodynamicsConfListener;
import edu.duke.cs.osprey.design.models.MoleculeDesign;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.BigMath;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

@Parameters(commandDescription = CommandPartitionFunction.CommandDescription)
public class CommandPartitionFunction extends RunDesignCommand {

    public static final String CommandName = "stability";
    static final String CommandDescription = "Estimate the partition function value(s) of different conformations";
    private final List<CommandAnalysis> confListeners = new LinkedList<>();

    @Parameter(names = "--thermodynamics", description = "Calculate the enthalpy and entropy of ensembles.")
    private boolean captureThermodynamics;

    @SuppressWarnings("MismatchedQueryAndUpdateOfCollection") // Ignored because updated by command line arguments
    @Parameter(names = "--energy", description = "Analyze the energy of conformation(s).")
    private List<Integer> captureEnergies = new ArrayList<>();

    @Parameter(names = "--max-num-confs", description = "Sets an upper bound on the number of conformations evaluated.")
    private int maxNumberConfs = -1;

    @Parameter(names = "--design-info", description = "Print information about the design and exit")
    private boolean printDesignInfo;

    private ConfEnergyCalculator confEnergyCalc;
    private PartitionFunction pFunc;
    private RCs rcs;

    @Override
    public int run(JCommander commander, String[] args) {

        var retVal = processHelpAndNoArgs(commander, args);

        if (retVal.isPresent()) {
            return retVal.get();
        }


        var designOpt = parseDesignSpec(MoleculeDesign.class);
        if (designOpt.isEmpty()) {
            return Main.Failure;
        }

        if (printDesignInfo) {
            return printDesignDebugInfo(designOpt.get());
        }

        return runStabilityDesign(designOpt.get());
    }

    @Override
    public String getCommandName() {
        return CommandName;
    }

    @Override
    public String getCommandDescription() {
        return CommandDescription;
    }

    private int printDesignDebugInfo(MoleculeDesign design) {
        var confSpace = delegate.createConfSpace(design.molecule, new ForcefieldParams());
        var numConfs = confSpace.getNumConformations();
        System.out.printf("Design: %s%n", design.designName);
        System.out.printf("Epsilon: %f%n", delegate.epsilon);
        System.out.printf("Number of conformations in design:\t%s%n", numConfs.toString());
        return Main.Success;
    }

    private int runStabilityDesign(MoleculeDesign design) {
        /* This reads parm96a.dat, which contains the energy parameters of DNA, RNA, and protein residues */
        var ffParams = new ForcefieldParams();

        /* Maintains flexibility information with the molecule, and can use that to make new molecules */
        var confSpace = delegate.createConfSpace(design.molecule, ffParams);

        if (delegate.verifyInput) {
            return Main.Success;
        }

        /* Decides whether to use CPU(s) and/or GPU(s) (purely implementation specific) */
        var parallelism = delegate.getParallelism();

        /* Used to calculate energies of a molecule, also used to minimize the molecule */
        try (var energyCalculator = new EnergyCalculator.Builder(confSpace, ffParams)
                .setParallelism(parallelism).build()) {

            confEnergyCalc = new ConfEnergyCalculator.Builder(confSpace, energyCalculator)
                    .build();

            /* Contains the confSpace and a pruning matrix */
            rcs = new RCs(confSpace);

            var epsilon = delegate.epsilon > 0 ? delegate.epsilon : 0.63;

            var energyMatrix = new SimplerEnergyMatrixCalculator.Builder(confEnergyCalc)
                    .build()
                    .calcEnergyMatrix();
            var lowerAStarTree = new ConfAStarTree.Builder(energyMatrix, rcs)
                    .setMPLP()
                    .build();
            var upperAStarTree = new ConfAStarTree.Builder(energyMatrix, rcs)
                    .setMPLP()
                    .build();

            try (var ctx = energyCalculator.tasks.contextGroup()) {
                pFunc = new GradientDescentPfunc(confEnergyCalc, lowerAStarTree, upperAStarTree, rcs.getNumConformations());
                pFunc.init(epsilon);

                pFunc.setInstanceId(0);
                pFunc.putTaskContexts(ctx);

                addListeners();
                pFunc.compute(maxNumberConfs > 0 ? maxNumberConfs : Integer.MAX_VALUE);
            }

            printResults();
            return Main.Success;
        }
    }

    private void printResults() {
        var numberFormat = NumberFormat.getPercentInstance();
        var percentEvaluated = numberFormat.format(new BigMath(PartitionFunction.decimalPrecision).set(pFunc.getNumConfsEvaluated()).div(rcs.getNumConformations().doubleValue()).get());

        System.out.println(String.format("Evaluated %s of conf space (%d / %s)", percentEvaluated, pFunc.getNumConfsEvaluated(), rcs.getNumConformations().toString()));
        System.out.println(pFunc.makeResult());

        for (var listener : confListeners) {
            listener.printResults();
        }
    }

    private void addListeners() {
        if (!captureEnergies.isEmpty()) {
            final var oneIndexed = captureEnergies.stream().map(x -> x - 1).collect(Collectors.toList());
            final var listener = new EnergyAnalysisConfListener(confEnergyCalc, oneIndexed);
            confListeners.add(listener);
            pFunc.setConfListener(listener);
        }

        if (captureThermodynamics) {
            final var listener = new ThermodynamicsConfListener();
            confListeners.add(listener);
            pFunc.setConfListener(listener);
        }
    }
}
