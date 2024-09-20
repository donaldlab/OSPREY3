package edu.duke.cs.osprey.design.commands;

// from Nate (change Kstar to *)
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.io.ScoredSequenceCsvPrinter;
import edu.duke.cs.osprey.design.io.ScoredSequenceEnsembleWriter;
import edu.duke.cs.osprey.design.io.ScoredSequencePFuncStatsWriter;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.ematrix.compiled.ErefCalculator;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.*;
import edu.duke.cs.osprey.kstar.pfunc.NewGradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

// from BBKstar
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import org.eclipse.collections.impl.block.function.MaxSizeFunction;

@Parameters(commandDescription = CompiledConfSpaceBBKStar.CommandDescription)
public class CompiledConfSpaceBBKStar implements CliCommand {

    public static final String CommandName = "bbkstar";
    public static final String CommandDescription = "Run the Branch and Bound K* algorithm with compiled conformation spaces.";

    @Parameter(names = "--complex-confspace", description = "Path to the compiled complex conformation space file.", required = true)
    private String complexConfSpacePath;

    @Parameter(names = "--design-confspace", description = "Path to the compiled design conformation space file.", required = true)
    private String designConfSpacePath;

    @Parameter(names = "--target-confspace", description = "Path to the compiled target conformation space file.", required = true)
    private String targetConfSpacePath;

    @Parameter(names = "--stability-threshold", description = "Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence. Set to a negative number to disable.")
    public double stabilityThreshold = 5.0;

    @Parameter(names = "--ensemble-dir", description = "Directory to save each sequence's structural ensemble in. Defaults to current working directory.")
    public String ensembleDir = ".";

    @Parameter(names = "--write-n-confs", description = "The number (n) of best conformations to write in each sequence's ensemble. Defaults to 10.")
    public int writeNConfs = 10;

    @Parameter(names = "--max-simultaneous-mutations", description = "How many residues should concurrently be allowed to mutate to non-wildtype amino acids?")
    public int maxSimultaneousMutations = 2;

    @Parameter(description = "The approximation accuracy. Z* = (1 - epsilon)Z. Values closer to 0 improve approximation accuracy.", names={"--epsilon", "-e"})
    double epsilon = 0.683;

    @Parameter(names = "--num-best-seqs", description = "The number of sequences before BBKStar stops.")
    public int NumBestSequences = 20;

    @Parameter(names = "--show-pfunc-prog", description = "Output the partition function progress during search?")
    boolean showPfunc = false;

    @Parameter(names = "--num-confs-batch", description = "Number of conformations per batch. Make this an even multiple of available threads.")
    public int numConfsBatch = 8;

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

        var start = System.currentTimeMillis();

        var design = ConfSpace.fromBytes(FileTools.readFileBytes(designConfSpacePath));
        var complex = ConfSpace.fromBytes(FileTools.readFileBytes(complexConfSpacePath));
        var target = ConfSpace.fromBytes(FileTools.readFileBytes(targetConfSpacePath));

        var taskExecutor = new Parallelism(Runtime.getRuntime().availableProcessors(), 0, 0)
                .makeTaskExecutor();

        var end1 = System.currentTimeMillis();
        System.out.printf("Took %f seconds to get files and start task%n", (end1 - start) / 1000.);

        // format Kstar score information
        KStarScoreWriter.Formatter testFormatter = info ->
                String.format("%3d %s   protein: %s   ligand: %s   complex: %s   K*: %s",
                        info.sequenceNumber,
                        info.sequence.toString(Sequence.Renderer.ResType),
                        info.kstarScore.protein.toString(),
                        info.kstarScore.ligand.toString(),
                        info.kstarScore.complex.toString(),
                        info.kstarScore.toString()
                );

        // customize Kstar and BBKstar settings
        KStarSettings kstarSettings = new KStarSettings.Builder()
                .setEpsilon(epsilon)
                .setStabilityThreshold(stabilityThreshold)
                // change this value for large designs
                .setMaxSimultaneousMutations(maxSimultaneousMutations)
                // set false to disable lots of printouts
                .setShowPfuncProgress(showPfunc)
                .addScoreConsoleWriter(testFormatter)
                .build();
        BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
                 // # of best seqs before BBK* stops
                 .setNumBestSequences(NumBestSequences)
                // make this even multiple of available threads
                .setNumConfsPerBatch(numConfsBatch)
                .build();
        BBKStar bbkstar = new BBKStar(target, design, complex, kstarSettings, bbkstarSettings);


        var end2 = System.currentTimeMillis();
        System.out.printf("Took %f seconds to do Kstar/BBKstar settings%n", (end2 - end1) / 1000.);

        for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

            // pass ConfSpace info
            ConfSpace confSpace = (ConfSpace) info.confSpace;

            // determine how residue interactions distributed among fragments
            PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
            boolean includeStaticStatic = true;
            ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(confSpace);

            // calculate minimized reference energies
            SimpleReferenceEnergies eref = new ErefCalculator.Builder(ecalc)
                    .setMinimize(true)
                    .build()
                    .calc();

            // create energy calculator amenable to ccsx file format
            info.confEcalcMinimized = new ConfEnergyCalculatorAdapter.Builder(ecalc, taskExecutor)
                    .setPosInterDist(posInterDist)
                    .setReferenceEnergies(eref)
                    .setMinimize(true)
                    .setIncludeStaticStatic(includeStaticStatic)
                    .build();

            // build energy matrix + A* search tree
            EnergyMatrix ematMinimized = new EmatCalculator.Builder(ecalc)
                    .setPosInterDist(posInterDist)
                    .setReferenceEnergies(eref)
                    .setMinimize(true)
                    .setIncludeStaticStatic(includeStaticStatic)
                    .build()
                    .calc();
            info.confSearchFactoryMinimized = (rcs) ->
                    new ConfAStarTree.Builder(ematMinimized, rcs)
                            .setTraditional()
                            .build();

            // BBK* needs rigid energies too
            EnergyMatrix ematRigid = new EmatCalculator.Builder(ecalc)
                    .setPosInterDist(posInterDist)
                    .setReferenceEnergies(eref)
                    .setMinimize(false)
                    .setIncludeStaticStatic(includeStaticStatic)
                    .build()
                    .calc();
            info.confSearchFactoryRigid = (rcs) ->
                    new ConfAStarTree.Builder(ematRigid, rcs)
                            .setTraditional()
                            .build();

            // use gradient descent to min pfuncs
            info.pfuncFactory = rcs -> new GradientDescentPfunc(
                    info.confEcalcMinimized,
                    info.confSearchFactoryMinimized.make(rcs),
                    info.confSearchFactoryMinimized.make(rcs),
                    rcs.getNumConformations()
            ).setPreciseBcalc(true);
        }

        var end3 = System.currentTimeMillis();
        System.out.printf("Took %f seconds to define confspace%n", (end3 - end2) / 1000.);

        // run BBK*
        List<ScoredSequence> sequences = bbkstar.run(taskExecutor);

        var end4 = System.currentTimeMillis();
        System.out.printf("Took %f seconds to run BBKStar%n", (end4 - end3) / 1000.);

        int lsize = sequences.size();
        System.out.println("Length of sequences: " + lsize);

        // make Sequence Analyzer + save PDB ensembles
        SequenceAnalyzer analyzer = new SequenceAnalyzer(bbkstar);
        int counter = 1;
        System.out.println("BBKstar results: ");
        for (ScoredSequence sequence : sequences) {
            System.out.println(" " + counter + "," + sequence.sequence() + "," + sequence.score());
            counter++;

            // set # conformations printed in ensemble + analyze
            SequenceAnalyzer.Analysis analysis = analyzer.analyze(sequence.sequence(), writeNConfs);

            // formats seqstr for file outputs (only changes filename)
            String seqstr = sequence.sequence().toString(Sequence.Renderer.ResTypeMutations)
                    .replace(' ', '-');

            // set where file is saved + name
            File ensembleFile = new File(String.format(ensembleDir + "/seq.%s.pdb", seqstr));

            // write the PDB
            analysis.writePdb(ensembleFile.getAbsolutePath(), String.format("Top %d conformations for sequence %s",
                    writeNConfs, sequence.sequence()));
        }

        var end5 = System.currentTimeMillis();
        System.out.printf("Took %f seconds for sequence analyzer%n", (end5 - end4) / 1000.);

        // cleanup
        for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {
            if (info.confEcalcMinimized != null) {
                ((ConfEnergyCalculatorAdapter) info.confEcalcMinimized).confEcalc.close();
            }
        }

        var end6 = System.currentTimeMillis();
        System.out.printf("Took %f seconds to cleanup%n", (end6 - end5) / 1000.);



        return Main.Success;
    }
}
