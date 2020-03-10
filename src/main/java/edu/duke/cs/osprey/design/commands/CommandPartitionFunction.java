package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.analysis.CommandAnalysis;
import edu.duke.cs.osprey.design.analysis.EnergyAnalysisConfListener;
import edu.duke.cs.osprey.design.analysis.ThermodynamicsConfListener;
import edu.duke.cs.osprey.design.models.StabilityDesign;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.tools.BigMath;

import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

@Parameters(commandDescription = CommandPartitionFunction.CommandDescription)
public class CommandPartitionFunction extends RunnableCommand {

    public static final String CommandName = "stability";
    static final String CommandDescription = "Estimate the partition function value(s) of different conformations";
    private final List<CommandAnalysis> confListeners = new LinkedList<>();

    @Parameter(names = "--thermodynamics", description = "Calculate the enthalpy and entropy of ensembles.")
    private boolean captureThermodynamics;

    @SuppressWarnings("MismatchedQueryAndUpdateOfCollection") // Ignored because updated by command line arguments
    @Parameter(names = "--energy", description = "Analyze the energy of conformation(s).")
    private List<Long> captureEnergies = new ArrayList<>();

    @Parameter(names = "--max-num-confs", description = "Sets an upper bound on the number of conformations evaluated.")
    private int maxNumberConfs = -1;

    @Parameter(names = "--design-info", description = "Print information about the design and exit")
    private boolean printDesignInfo;

    private ConfEnergyCalculator confEnergyCalc;

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

        if (printDesignInfo) {
            return printDesignDebugInfo(design);
        }

        return runStabilityDesign(design);
    }

    @Override
    public String getCommandName() {
        return CommandName;
    }

    @Override
    public String getCommandDescription() {
        return CommandDescription;
    }

    private int printDesignDebugInfo(StabilityDesign design) {
        var confSpace = delegate.createConfSpace(design.molecule, new ForcefieldParams());
        var numConfs = confSpace.getNumConformations();
        System.out.println(String.format("Design: %s", design.designName));
        System.out.println(String.format("Epsilon: %f", design.epsilon));
        System.out.println(String.format("Number of conformations in design:\t%s", numConfs.toString()));
        return Main.Success;
    }

    private int runStabilityDesign(StabilityDesign design) {
        /* This reads parm96a.dat, which contains the energy parameters of DNA, RNA, and protein residues */
        var ffParams = new ForcefieldParams();

        /* Maintains flexibility information with the molecule, and can use that to make new molecules */
        var confSpace = delegate.createConfSpace(design.molecule, ffParams);

        if (delegate.verifyInput) {
            System.out.println("The stability design input file is valid.");
            return Main.Success;
        }

        /* Decides whether to use CPU(s) and/or GPU(s) (purely implementation specific) */
        var parallelism = delegate.getParallelism();

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
        // https://github.com/donaldlab/OSPREY3/blob/sharkstar/test/edu/duke/cs/osprey/sharkstar/TestSHARKStarBound.java#L65
        confEnergyCalc = new ConfEnergyCalculator.Builder(confSpace, energyCalculator)
                .build();

        final var sequences = confSpace.seqSpace.getSequences();
        final var epsilon = delegate.epsilon > 0 ? delegate.epsilon : design.epsilon;

        for (Sequence seq : sequences) {
            var rc = seq.makeRCs(confSpace);
            var partitionFnBuilder = new PartitionFunctionFactory(confSpace, confEnergyCalc, design.designName);
            partitionFnBuilder.setUseGradientDescent();
            var pfunc = partitionFnBuilder.makePartitionFunctionFor(rc, epsilon);
            addListeners(pfunc, seq.toString());
            pfunc.compute(maxNumberConfs > 0 ? maxNumberConfs : Integer.MAX_VALUE);
            printResults(seq, rc, pfunc);
        }

        return Main.Success;
    }

    private void printResults(Sequence seq, RCs rc, PartitionFunction pf) {
        var numberFormat = NumberFormat.getPercentInstance();
        var percentEvaluated = numberFormat.format(new BigMath(PartitionFunction.decimalPrecision).set(pf.getNumConfsEvaluated()).div(rc.getNumConformations().doubleValue()).get());
        System.out.println(String.format("Evaluated %s of conf space (%d / %s)", percentEvaluated, pf.getNumConfsEvaluated(), rc.getNumConformations().toString()));
        System.out.println(seq.toString(Sequence.Renderer.AssignmentMutations));
        System.out.println(pf.makeResult());

        for (var listener : confListeners) {
            listener.printResults();
        }
    }

    private void addListeners(PartitionFunction pfunc, String sequenceDescription) {

        for (var listener : delegate.makeListeners(pfunc, confEnergyCalc, sequenceDescription)) {
            pfunc.addConfListener(listener);
            confListeners.add(listener);
        }

        if (!captureEnergies.isEmpty()) {
            final var oneIndexed = captureEnergies.stream().map(x -> x - 1).collect(Collectors.toList());
            final var listener = new EnergyAnalysisConfListener(confEnergyCalc, oneIndexed);
            confListeners.add(listener);
            pfunc.addConfListener(listener);
        }

        if (captureThermodynamics) {
            final var listener = new ThermodynamicsConfListener();
            confListeners.add(listener);
            pfunc.addConfListener(listener);
        }
    }
}

