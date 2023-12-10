package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.*;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.io.ScoredSequenceCsvPrinter;
import edu.duke.cs.osprey.design.io.ScoredSequenceEnsembleWriter;
import edu.duke.cs.osprey.design.io.ScoredSequencePFuncStatsWriter;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.ematrix.compiled.ErefCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.gpu.Precision;
import edu.duke.cs.osprey.kstar.KStarSettings;
import edu.duke.cs.osprey.kstar.NewKStar;
import edu.duke.cs.osprey.kstar.ScoredSequence;
import edu.duke.cs.osprey.kstar.pfunc.NewGradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.PrintWriter;
import java.util.List;

@Parameters(commandDescription = CompiledConfSpaceKStar.CommandDescription)
public class CompiledConfSpaceKStar implements CliCommand {

    public static final String CommandName = "kstar";
    public static final String CommandDescription = "Run the K* algorithm with compiled conformation spaces.";

    @Parameter(names = "--complex-confspace", description = "Path to the compiled complex conformation space file.", required = true)
    private String complexConfSpacePath;

    @Parameter(names = "--design-confspace", description = "Path to the compiled design conformation space file.", required = true)
    private String designConfSpacePath;

    @Parameter(names = "--target-confspace", description = "Path to the compiled target conformation space file.", required = true)
    private String targetConfSpacePath;

    @Parameter(names = "--stability-threshold", description = "Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence. Set to a negative number to disable.")
    public double stabilityThreshold = 5.0;

    @Parameter(names = {"--max-confs"}, description = "Number of lowest energy conformations to save. Off by default.")
    public int maxConfs = -1;

    @Parameter(names = "--ensemble-dir", description = "Directory to save each sequence's structural ensemble in. Defaults to current working directory.")
    public String ensembleDir = ".";

    @Parameter(names = "--write-n-confs", description = "The number (n) of best conformations to write in each sequence's ensemble. Defaults to 10.")
    public int writeNConfs = 10;

    @Parameter(names = "--max-simultaneous-mutations", description = "How many residues should concurrently be allowed to mutate to non-wildtype amino acids?")
    public int maxSimultaneousMutations = 2;

    @Parameter(description = "The approximation accuracy. Z* = (1 - epsilon)Z. Values closer to 0 improve approximation accuracy.", names={"--epsilon", "-e"})
    double epsilon = 0.683;

    @Parameter(names = "--precision", description = "The precision of the energy calculator's floating points { fp32 or fp64 }.", converter = PrecisionConverter.class)
    Precision precision = Precision.Float64;

    @Override
    public String getCommandName() {
        return CommandName;
    }

    @Override
    public String getCommandDescription() {
        return CommandDescription;
    }

    @Override
    public int run(JCommander commander, String[] args) {

        var design = ConfSpace.fromBytes(FileTools.readFileBytes(designConfSpacePath));
        var designConfCalc = ConfEnergyCalculator.makeBest(design, precision);

        // use a CPU task executor only for calculating reference energies
        var taskExecutor = new Parallelism(Runtime.getRuntime().availableProcessors(), 0, 0)
                .makeTaskExecutor();
        var designRefEnergies = new ErefCalculator.Builder(designConfCalc)
                .build()
                .calc(taskExecutor);

        System.out.println("Design ref energies");
        System.out.println(designRefEnergies.toString(design));

        var complex = ConfSpace.fromBytes(FileTools.readFileBytes(complexConfSpacePath));
        var complexConfCalc = ConfEnergyCalculator.makeBest(complex, precision);
        var complexRefEnergies = new ErefCalculator.Builder(complexConfCalc)
                .build()
                .calc(taskExecutor);

        System.out.println("Complex ref energies");
        System.out.println(complexRefEnergies.toString(complex));

        var target = ConfSpace.fromBytes(FileTools.readFileBytes(targetConfSpacePath));
        var targetConfCalc = ConfEnergyCalculator.makeBest(target, precision);
        var targetRefEnergies = new ErefCalculator.Builder(targetConfCalc)
                .build()
                .calc(taskExecutor);
        System.out.println("Target ref energies");
        System.out.println(targetRefEnergies.toString(target));

        var complexEnergyMatrix = new EmatCalculator.Builder(complexConfCalc)
                .setReferenceEnergies(complexRefEnergies)
                .build()
                .calc(taskExecutor);
        var targetEnergyMatrix = new EmatCalculator.Builder(targetConfCalc)
                .setReferenceEnergies(targetRefEnergies)
                .build()
                .calc(taskExecutor);
        var designEnergyMatrix = new EmatCalculator.Builder(designConfCalc)
                .setReferenceEnergies(designRefEnergies)
                .build()
                .calc(taskExecutor);

        var complexPosInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, complexRefEnergies);
        var targetPosInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, targetRefEnergies);
        var designPosInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, designRefEnergies);

        var settings = new KStarSettings.Builder()
                .setStabilityThreshold(stabilityThreshold)
                .setMaxNumConf(maxConfs > 0 ? maxConfs : Integer.MAX_VALUE)
                .setMaxSimultaneousMutations(maxSimultaneousMutations)
                .setEpsilon(epsilon)
                .build();

        var kstar = new NewKStar(target, design, complex, settings);
        kstar.protein.pfuncFactory = (rcs) -> new NewGradientDescentPfunc(
                targetConfCalc,
                new ConfAStarTree.Builder(targetEnergyMatrix, rcs).setTraditional().build(),
                new ConfAStarTree.Builder(targetEnergyMatrix, rcs).setTraditional().build(),
                rcs.getNumConformations(),
                targetPosInterGen,
                taskExecutor
        );
        kstar.protein.confEcalc = targetConfCalc;

        kstar.ligand.pfuncFactory = (rcs) -> new NewGradientDescentPfunc(
                designConfCalc,
                new ConfAStarTree.Builder(designEnergyMatrix, rcs).setTraditional().build(),
                new ConfAStarTree.Builder(designEnergyMatrix, rcs).setTraditional().build(),
                rcs.getNumConformations(),
                designPosInterGen,
                taskExecutor
        );
        kstar.ligand.confEcalc = designConfCalc;

        kstar.complex.pfuncFactory = (rcs) -> new NewGradientDescentPfunc(
                complexConfCalc,
                new ConfAStarTree.Builder(complexEnergyMatrix, rcs).setTraditional().build(),
                new ConfAStarTree.Builder(complexEnergyMatrix, rcs).setTraditional().build(),
                rcs.getNumConformations(),
                complexPosInterGen,
                taskExecutor
        );
        kstar.complex.confEcalc = complexConfCalc;

        var printWriter = new PrintWriter(System.out);
        var csvPrinter = new ScoredSequenceCsvPrinter(printWriter);
        var statsWriter = new ScoredSequencePFuncStatsWriter(ensembleDir);
        var ensembleWriter = new ScoredSequenceEnsembleWriter(writeNConfs, ensembleDir, complexRefEnergies, complexConfCalc, complex);

        kstar.putSequenceComputedListener(csvPrinter);
        kstar.putSequenceComputedListener(ensembleWriter);
        kstar.putSequenceComputedListener(statsWriter);

        List<ScoredSequence> sequences = kstar.run(taskExecutor);
        return Main.Success;
    }

    public class ValidatePrecision implements IParameterValidator {
        @Override
        public void validate(String name, String value) throws ParameterException {
            switch (value) {
                case "fp32":
                case "fp64":
                    break;
                default:
                    throw new ParameterException(String.format("Parameter %s must be fp32 or fp64", name));

            }
        }
    }

    public class PrecisionConverter implements IStringConverter<Precision> {
        @Override
        public Precision convert(String s) {
            return switch (s) {
                case "fp64" -> Precision.Float64;
                case "fp32" -> Precision.Float32;
                default -> throw new ParameterException(String.format("Parameter %s must be fp32 or fp64", "precision"));
            };
        }
    }
}
