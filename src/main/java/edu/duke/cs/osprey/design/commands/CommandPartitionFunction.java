package edu.duke.cs.osprey.design.commands;


import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.analysis.CommandAnalysis;
import edu.duke.cs.osprey.design.models.AminoAcid;
import edu.duke.cs.osprey.design.models.StabilityDesign;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.design.analysis.EnergyAnalysisConfListener;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.design.analysis.ThermodynamicsConfListener;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.BigMath;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaDeviceProp;

import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

import static jcuda.runtime.JCuda.cudaGetDeviceCount;
import static jcuda.runtime.JCuda.cudaGetDeviceProperties;

@Parameters(commandDescription = CommandPartitionFunction.CommandDescription)
public class CommandPartitionFunction extends RunnableCommand {

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

    @Parameter(names = {"--cuda", "-c"})
    private boolean useCuda;

    @Parameter(names = "--verify-design", description = "Verifies input parameters, but does not run the design")
    private boolean verifyInput;

    private ConfEnergyCalculator confEnergyCalc;
    private PartitionFunction pFunc;
    private RCs rcs;
    /* This reads parm96a.dat, which contains the energy parameters of DNA, RNA, and protein residues */
    private ForcefieldParams ffParams = new ForcefieldParams();

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

    private SimpleConfSpace createConfSpace(StabilityDesign design) {
        /* Reads a PDB file into a Molecule. */
        var molecule = PDBIO.read(design.molecule);

        /* Reads the templates, rotamers, and entropy info for a given forcefield */
        /*
            "/config/parm96a.dat",
            "/config/all_amino94.in",
            "/config/all_aminont94.in",
            "/config/all_aminoct94.in",
            "/config/all_nuc94_and_gr.in",
         */
        var extraTemplates = design.extraTemplates;
        var extraTemplateCoords = design.extraTemplatesCoordinates;

        var templateLibrary = new ResidueTemplateLibrary.Builder(ffParams.forcefld)
                .addTemplates(extraTemplates)
                .addTemplateCoords(extraTemplateCoords)
                .build();

        /* Strands combine a Molecule with design flexibility and templates */
        var strandBuilder = new Strand.Builder(molecule)
                .setTemplateLibrary(templateLibrary);

        /* Add in flexibility and mutability parameters */
        for (var mod : design.residueModifiers) {
            var identifier = mod.identity.positionIdentifier();
            strandBuilder.setResidueMutability(
                    identifier,
                    mod.mutability.stream().map(AminoAcid::toValue).collect(Collectors.toUnmodifiableList()),
                    mod.flexibility.includeStructureRotamer,
                    mod.flexibility.useContinuous
            );
        }
        var protein = strandBuilder.build();

        // Validate that design intentions match input structure
        var identities = design.residueModifiers.stream().map(m -> m.identity).collect(Collectors.toUnmodifiableList());
        for(var identity : identities) {
            var flex = protein.flexibility.get(identity.positionIdentifier());
            if (!flex.wildType.substring(0, 2).equalsIgnoreCase(identity.aminoAcidType.toValue().substring(0, 2))) {
                throw new RuntimeException(String.format("Design parameter thought residue %s was %s, but in structure is %s", identity.positionIdentifier(), identity.aminoAcidType.toValue(), flex.wildType));
            }
        }

        if (verifyInput) {
            return Main.Success;
        }

        /* Maintains flexibility information with the molecule, and can use that to make new molecules */

        return new SimpleConfSpace.Builder()
                .addStrand(protein)
                .setShellDistance(0)
                .build();
    }

    private int printDesignDebugInfo(StabilityDesign design) {
        var confSpace = createConfSpace(design);
        var numConfs = confSpace.getNumConformations();
        System.out.println(String.format("Design: %s", design.designName));
        System.out.println(String.format("Epsilon: %f", design.epsilon));
        System.out.println(String.format("Number of conformations in design:\t%s", numConfs.toString()));
        return Main.Success;
    }

    private int runStabilityDesign(StabilityDesign design) {
        /* Maintains flexibility information with the molecule, and can use that to make new molecules */
        var confSpace = createConfSpace(design);

        /* Decides whether to use CPU(s) and/or GPU(s) (purely implementation specific) */
        var parallelism = getParallelism();

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

        /* Contains the confSpace and a pruning matrix */
        rcs = new RCs(confSpace);

        var partitionFnBuilder = new PartitionFunctionFactory(confSpace, confEnergyCalc, "default");
        partitionFnBuilder.setUseGradientDescent();
        pFunc = partitionFnBuilder.makePartitionFunctionFor(rcs, delegate.epsilon > 0 ? delegate.epsilon : design.epsilon);
        addListeners();
        pFunc.compute(maxNumberConfs > 0 ? maxNumberConfs : Integer.MAX_VALUE);

        printResults();
        return Main.Success;
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
            pFunc.addConfListener(listener);
        }

        if (captureThermodynamics) {
            final var listener = new ThermodynamicsConfListener();
            confListeners.add(listener);
            pFunc.addConfListener(listener);
        }
    }

    private int getNumGpu() {
        var deviceCount = new int[]{0};
        cudaGetDeviceCount(deviceCount);
        return deviceCount[0];
    }

    private int getNumberGpuMultiprocessors(int deviceIdx) {
        var deviceProps = new cudaDeviceProp();
        cudaGetDeviceProperties(deviceProps, deviceIdx);
        return deviceProps.multiProcessorCount;
    }

    private Parallelism getParallelism() {
        if (useCuda) {
            JCuda.setExceptionsEnabled(true);
            var numGpus = getNumGpu();
            var numMultiprocessors = getNumberGpuMultiprocessors(0);
            return new Parallelism(Runtime.getRuntime().availableProcessors(), numGpus, numMultiprocessors);
        }

        return new Parallelism(Runtime.getRuntime().availableProcessors(), 0, 0);
    }
}
