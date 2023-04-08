package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.google.common.collect.Lists;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.ematrix.compiled.ErefCalculator;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.*;
import edu.duke.cs.osprey.kstar.pfunc.NewGradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.File;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

@Parameters(commandDescription = CommandCKStar.CommandDescription)
public class CommandCKStar extends RunnableCommand {

    public static final String CommandName = "ckstar";
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

        var start = System.currentTimeMillis();

//        cleanupStuffFromPreviousRuns();

        var numConfs = 200;

        var complex = ConfSpace.fromBytes(FileTools.readFileBytes(complexConfSpacePath));
        var design = ConfSpace.fromBytes(FileTools.readFileBytes(designConfSpacePath));
        var target = ConfSpace.fromBytes(FileTools.readFileBytes(targetConfSpacePath));
        var parallelism = delegate.getParallelism();
//        var parallelism = new Parallelism(16, 0, 0);
        var taskExecutor = parallelism.makeTaskExecutor();

        var settings = new NewKStar.Settings.Builder()
                .setStabilityThreshold(stabilityThreshold)
                .setExternalMemory(false)
                .setMaxNumConf(numConfs)
                .resume(true)
                .addScoreFileWriter(Paths.get("sequences.tsv").toFile())
                .build();

        var kstar = new NewKStar(target, design, complex, settings);

        for (var confSpaceInfo : kstar.confSpaceInfos()) {
            var energyCalculator = new CPUConfEnergyCalculator((ConfSpace) confSpaceInfo.confSpace);
            confSpaceInfo.confEcalc = energyCalculator;

            var referenceEnergies = new ErefCalculator.Builder(energyCalculator)
                    .build()
                    .calc(taskExecutor);

            var energyMatrix = new EmatCalculator.Builder(energyCalculator)
                    .setReferenceEnergies(referenceEnergies)
                    .setCacheFile(new File(String.format("emat.%s.dat", confSpaceInfo.id)))
                    .build()
                    .calc(taskExecutor);

            PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, referenceEnergies);

            confSpaceInfo.pfuncFactory =
                    (rcs, ctxGroup) -> new NewGradientDescentPfunc(
                            energyCalculator,
                            new ConfAStarTree.Builder(energyMatrix, rcs).setTraditional().build(),
                            new ConfAStarTree.Builder(energyMatrix, rcs).setTraditional().build(),
                            rcs.getNumConformations(),
                            posInterGen,
                            taskExecutor,
                            ctxGroup
                    );
        }

        List<NewKStar.ScoredSequence> sequences = kstar.run(taskExecutor);

        int i = 0;
        for (var sequence : sequences) {
            System.out.println(sequence);
        }

        var stop = System.currentTimeMillis();
        System.out.printf("Took %f seconds to run%n", (stop - start) / 1000.);

        /*
        SequenceAnalyzer analyzer = new SequenceAnalyzer(Lists.newArrayList(kstar.confSpaceInfos()));

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
