package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.directors.KStarDirector;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.tools.FileTools;

import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

@Parameters(commandDescription = CommandKStarBoundedMem.CommandDescription)
public class CommandKStarBoundedMem extends RunnableCommand {

    public static final String CommandName = "kstar-bounded-mem";
    public static final String CommandDescription = "Run K* using the bounded-memory energy calculator";

    @Parameter(names = "--complex-confspace", description = "Path to the compiled complex conformation space file.", required = true)
    private String complexConfSpacePath;

    @Parameter(names = "--design-confspace", description = "Path to the compiled design conformation space file.", required = true)
    private String designConfSpacePath;

    @Parameter(names = "--target-confspace", description = "Path to the compiled target conformation space file.", required = true)
    private String targetConfSpacePath;

    @Parameter(names = "--stability-threshold", description = "Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence. Set to a negative number to disable.")
    public double stabilityThreshold = 5.0;

    @Parameter(names = "--g-width-max", description = "Set the precision for calculated free energy values. A small width is more precise, but can take longer to calculate.")
    public double gWidthMax = 0.1;
    @Parameter(names = "--report-state-progress", description = "Show progress for each free energy calculation for each sequence for each state.")
    public boolean reportStateProgress = true;

    @Parameter(names = "--track-ensembles", description = "Periodically write out a pdb ensemble of the 5 lowest energy structure for sequence to a directory.")
    public boolean trackEnsembles = true;

    @Parameter(names = "--track-ensembles-directory", description = "directory to write pdb ensembles to.")
    public String trackEnsemblesDirectory = ".";

    @Parameter(names = "--free-energy-mem-gib", description = "Gigabytes of memory to use for the free energy calculation.")
    public int freeEnergyGib = 16;

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
        var complex = ConfSpace.fromBytes(FileTools.readFileBytes(complexConfSpacePath));
        var design = ConfSpace.fromBytes(FileTools.readFileBytes(designConfSpacePath));
        var target = ConfSpace.fromBytes(FileTools.readFileBytes(targetConfSpacePath));
        var parallelism = delegate.getParallelism();

        var builder = new KStarDirector.Builder(complex, design, target);
        builder.setGWidthMax(gWidthMax);
        builder.setMaxSimultaneousMutations(delegate.maxSimultaneousMutations);
        builder.setStabilityThreshold(stabilityThreshold);
        builder.setReportStateProgress(reportStateProgress);
        if (trackEnsembles) {
            builder.setEnsembleTracking(5, Paths.get(trackEnsemblesDirectory).toFile());
        }

        KStarDirector director = builder.build();

        Coffee.Builder coffeeBuilder = new Coffee.Builder(director.confSpace)
                .setParallelism(parallelism);

        for (var state : director.confSpace.states) {
            var config = new Coffee.StateConfig(state);
            config.posInterGen = new PosInterGen(PosInterDist.TighterBounds, null);
            coffeeBuilder.configState(config);
        }

        Coffee coffee = coffeeBuilder.build();
        coffee.run(director);

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
