package edu.duke.cs.osprey.design.commands;

// from Nate (change Kstar to *)
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.ematrix.compiled.ErefCalculator;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.kstar.*;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.File;
import java.util.List;

// from BBKstar
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;

@Parameters(commandDescription = CommandKStar.CommandDescription)
public class CommandBBKStar extends DelegatingCommand {

    public static final String CommandName = "bbkstar";
    public static final String CommandDescription = "Run Branch and Bound K* using compiled conformation spaces.";

    @Parameter(names = "--complex-confspace", description = "Path to the compiled complex conformation space file.", required = true)
    private String complexConfSpacePath;

    @Parameter(names = "--design-confspace", description = "Path to the compiled design conformation space file.", required = true)
    private String designConfSpacePath;

    @Parameter(names = "--target-confspace", description = "Path to the compiled target conformation space file.", required = true)
    private String targetConfSpacePath;

    @Parameter(names = "--stability-threshold", description = "Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence. Set to a negative number to disable.")
    public double stabilityThreshold = 5.0;

    @Override
    public int run(JCommander commander, String[] args) {
        cleanupStuffFromPreviousRuns();
        var start = System.currentTimeMillis();

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

        var complex = ConfSpace.fromBytes(FileTools.readFileBytes(complexConfSpacePath));
        var design = ConfSpace.fromBytes(FileTools.readFileBytes(designConfSpacePath));
        var target = ConfSpace.fromBytes(FileTools.readFileBytes(targetConfSpacePath));
        var parallelism = delegate.getParallelism();
        var taskExecutor = parallelism.makeTaskExecutor();

        // set Kstar settings
        KStarSettings kstarsettings = new KStarSettings.Builder()
                .setStabilityThreshold(stabilityThreshold)
                .setShowPfuncProgress(true)
                .setMaxNumConf(delegate.maxConfs)
                .build();

        // set BBKStar settings
        BBKStar.Settings bbkstarsettings = new BBKStar.Settings.Builder()
                // # of best seqs before BBK* stops
                .setNumBestSequences(20)
                // make this even multiple of available threads
                .setNumConfsPerBatch(8)
                .build();

        // build BBKStar
        BBKStar bbkstar = new BBKStar(target, design, complex, kstarsettings, bbkstarsettings);

        // make E calc and confspace

        for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

            ConfSpace confSpace = (ConfSpace) info.confSpace;

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


        // run BBK* search
        List<ScoredSequence> sequences = bbkstar.run(taskExecutor);
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
            SequenceAnalyzer.Analysis analysis = analyzer.analyze(sequence.sequence(), delegate.writeNConfs);
            System.out.println(analysis);

            // formats seqstr for file outputs (only changes filename)
            String seqstr = sequence.sequence().toString(Sequence.Renderer.ResTypeMutations)
                    .replace(' ', '-');

            // set where file is saved + name
            File ensembleFile = new File(delegate.ensembleDir, String.format("seq.%s.pdb", seqstr));

            // write the PDB (can also set filepath here)
            analysis.writePdb(ensembleFile.getAbsolutePath(), String.format("Top %d conformations for sequence %s",
                    delegate.writeNConfs, sequence.sequence()));
        }

        var end = System.currentTimeMillis();
        System.out.printf("Took %f seconds%n", (end - start) / 1000.);


        // cleanup
        for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {
            if (info.confEcalcMinimized != null) {
                ((ConfEnergyCalculatorAdapter)info.confEcalcMinimized).confEcalc.close();
            }
        }

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
