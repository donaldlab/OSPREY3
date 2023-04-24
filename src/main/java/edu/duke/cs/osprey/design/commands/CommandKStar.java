package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.ematrix.compiled.ErefCalculator;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.tools.FileTools;

import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

@Parameters(commandDescription = CommandKStar.CommandDescription)
public class CommandKStar extends RunnableCommand {

    public static final String CommandName = "kstar";
    public static final String CommandDescription = "Run K* using compiled conformation spaces.";

    @Parameter(names = "--complex-confspace", description = "Path to the compiled complex conformation space file.", required = true)
    private String complexConfSpacePath;

    @Parameter(names = "--design-confspace", description = "Path to the compiled design conformation space file.", required = true)
    private String designConfSpacePath;

    @Parameter(names = "--target-confspace", description = "Path to the compiled target conformation space file.", required = true)
    private String targetConfSpacePath;

    @Parameter(names = "--stability-threshold", description = "Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence. Set to a negative number to disable.")
    public double stabilityThreshold = 5.0;

    private String[] args;

    private Map<String, String> dbSettings = new HashMap<>();
    private static final String dbHostnameKey = "dbHostname";
    private static final String dbPortKey = "dbPort";
    private static final String dbNameKey = "dbName";
    private static final String dbUserNameKey = "dbUserName";
    private static final String dbPasswordKey = "dbPassword";

    @Override
    public int run(JCommander commander, String[] args) {
        this.args = args;

        cleanupStuffFromPreviousRuns();

        var start = System.currentTimeMillis();
        var numConfs = 200;

        var complex = ConfSpace.fromBytes(FileTools.readFileBytes(complexConfSpacePath));
        var design = ConfSpace.fromBytes(FileTools.readFileBytes(designConfSpacePath));
        var target = ConfSpace.fromBytes(FileTools.readFileBytes(targetConfSpacePath));
        var parallelism = delegate.getParallelism();
//        var parallelism = new Parallelism(16, 0, 0);
        var taskExecutor = parallelism.makeTaskExecutor();

        var settings = new KStar.Settings.Builder()
                .setStabilityThreshold(stabilityThreshold)
                .setMaxNumConf(numConfs)
                .build();

        KStar kstar = new KStar(target, design, complex, settings);

        for (KStar.ConfSpaceInfo confSpaceInfo : kstar.confSpaceInfos()) {
            var energyCalculator = new CPUConfEnergyCalculator((ConfSpace) confSpaceInfo.confSpace);
            var referenceEnergies = new ErefCalculator.Builder(energyCalculator)
                    .build()
                    .calc(taskExecutor);

            var energyMatrix = new EmatCalculator.Builder(energyCalculator)
                    .setReferenceEnergies(referenceEnergies)
//                    .setCacheFile(new File(String.format("emat.%s.dat", confSpaceInfo.id)))
                    .build()
                    .calc(taskExecutor);

            var confEnergyCalculator = new ConfEnergyCalculatorAdapter.Builder(energyCalculator, taskExecutor)
                    .setReferenceEnergies(referenceEnergies)
                    .build();

            confSpaceInfo.confEcalc = confEnergyCalculator;

            confSpaceInfo.pfuncFactory = (rcs) -> new GradientDescentPfunc(
                    confSpaceInfo.confEcalc,
                    new ConfAStarTree.Builder(energyMatrix, rcs).setTraditional().build(),
                    new ConfAStarTree.Builder(energyMatrix, rcs).setTraditional().build(),
                    rcs.getNumConformations()
            );
        }

        List<KStar.ScoredSequence> sequences = kstar.run(taskExecutor);
        for (var sequence : sequences) {
            System.out.println(sequence);
        }

        var end = System.currentTimeMillis();
        System.out.printf("Took %f seconds%n", (end - start) / 1000.);

        /*
        SequenceAnalyzer analyzer = new SequenceAnalyzer(kstar);
        var seq = 1;
        for (KStar.ScoredSequence sequence : sequences) {
            System.out.println("result:");
            System.out.println("\tsequence: " + sequence.sequence);
            System.out.println("\tscore: " + sequence.score);

            SequenceAnalyzer.Analysis analysis = analyzer.analyze(sequence.sequence, numConfs);
            System.out.println(analysis);

            analysis.writePdb(String.format("ensemble-%d.pdb", seq++), String.format("Top %d conformations for sequence %s", numConfs, sequence.sequence));
        }
         */

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
