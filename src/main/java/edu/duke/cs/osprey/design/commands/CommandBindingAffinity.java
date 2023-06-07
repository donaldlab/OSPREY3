package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.confspace.StrandFlex;
import edu.duke.cs.osprey.design.*;
import edu.duke.cs.osprey.design.models.AffinityDesign;
import edu.duke.cs.osprey.design.models.MoleculeDto;
import edu.duke.cs.osprey.design.models.ResidueModifier;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.amber.ForcefieldFileParser;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarSettings;
import edu.duke.cs.osprey.kstar.ScoredSequence;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import one.util.streamex.IntStreamEx;
import one.util.streamex.StreamEx;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Stream;

@Parameters(commandDescription = CommandBindingAffinity.CommandDescription)
public class CommandBindingAffinity extends DelegatingCommand {

    public static final String CommandName = "affinity";
    public static final String CommandDescription = "Compute an epsilon approximation to binding affinity (K*).";

    @Parameter(names = "--do-scan", description = "Runs a scan using the scan settings specified in the design.")
    public boolean doScan;

    @Parameter(names = "--scan-flex-distance", description = "Distance (in angstroms) around mutable residues in which residues are flexible.")
    public double scanFlexDistance = 2.0;

    @Parameter(names = "--scan-distance", description = "Distance (in angstroms) around a scan's target residue in which residues are made mutable.")
    public double scanDistance = 5.0;

    @Parameter(names = "--scan-output", description = "Specifies the output directory to save the scan designs in.")
    public String scanOutput;

    @Parameter(names = "--max-num-confs", description = "Sets an upper bound on the number of conformations evaluated.")
    private int maxNumberConfs = -1;

    @Parameter(names = "--use-markstar", description = "Use MARK* instead of Gradient Descent for the partition function calculation.")
    public boolean useMarkstar;

    @Parameter(names = "--stability-threshold", description = "Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence. Set to a negative number to disable.")
    public double stabilityThreshold = 5.0;

    private String[] args;

    private Map<String, String> dbSettings = new HashMap<>();
    private static final String dbHostnameKey = "dbHostname";
    private static final String dbPortKey = "dbPort";
    private static final String dbNameKey = "dbName";
    private static final String dbUserNameKey = "dbUserName";
    private static final String dbPasswordKey = "dbPassword";
    private static final List<String> allMutableAaTypes = List.of("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIE", "HID", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL");

    @Override
    public int run(JCommander commander, String[] args) {
        this.args = args;

        return processHelpAndNoArgs(commander, args)
                .orElseGet(() -> parseAndValidate(delegate.design)
                .map(this::runAffinityDesign)
                .orElse(Main.Failure));
    }

    private int runAffinityDesign(AffinityDesign design) {
        var paramsAndStrands = new ForceFieldParamsAndStrands(delegate, design);

        if (doScan && design.scanSettings != null) {
            return makeScanDesigns(design, paramsAndStrands.makeAllResidues());
        }

        // Exit early if just trying to validate input
        if (delegate.verifyInput) {
            System.out.println("Design file validated.");
            return Main.Success;
        }

        /* Used to calculate energies of a molecule, also used to minimize the molecule */
        var minimizingECalc = new EnergyCalculator.Builder(paramsAndStrands.complexConfSpace, paramsAndStrands.forcefieldParams)
                .setParallelism(delegate.getParallelism())
                .build();
        var rigidECalc = new EnergyCalculator.SharedBuilder(minimizingECalc)
                .setIsMinimizing(false)
                .build();

        var epsilon = delegate.epsilon > 0 ? delegate.epsilon : 0.999999;
        var kstar = new KStar(paramsAndStrands.proteinConfSpace, paramsAndStrands.ligandConfSpace, paramsAndStrands.complexConfSpace, makeKStarSettings(epsilon));

        for (var info : kstar.confSpaceInfos()) {
            var referenceEnergies = new SimpleReferenceEnergies.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc)
                    .build();
            var minimizingConfECalc = new ConfEnergyCalculator.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc)
                    .setReferenceEnergies(referenceEnergies)
                    .build();
            info.confEcalc = minimizingConfECalc;

            var minimizedEnergyMatrix = new SimplerEnergyMatrixCalculator.Builder(minimizingConfECalc)
                    .build()
                    .calcEnergyMatrix();

            if (useMarkstar) {
                var rigidConfECalc = new ConfEnergyCalculator(info.confEcalc, rigidECalc);
                var rigidEnergymatrix = new SimplerEnergyMatrixCalculator.Builder(rigidConfECalc)
                        .build()
                        .calcEnergyMatrix();

                info.pfuncFactory = (rcs) -> {
                    var pfn = new MARKStarBoundFastQueues(minimizingConfECalc.confSpace, rigidEnergymatrix, minimizedEnergyMatrix, info.confEcalc, rcs, minimizingECalc.parallelism);
                    pfn.setCorrections(new UpdatingEnergyMatrix(info.confEcalc.confSpace, minimizedEnergyMatrix, info.confEcalc));
                    return pfn;
                };
            } else {
                info.pfuncFactory = (rcs) -> new GradientDescentPfunc(
                        info.confEcalc,
                        new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                        new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                        rcs.getNumConformations()
                );
            }
        }

        printResults(kstar.run(minimizingECalc.tasks));
        return Main.Success;
    }

    // I don't like this, think of a better way to encapsulate this
    public static class ForceFieldParamsAndStrands {
        private final SimpleConfSpace complexConfSpace;
        public ForcefieldParams forcefieldParams;
        public SimpleConfSpace proteinConfSpace;
        public SimpleConfSpace ligandConfSpace;
        public List<Strand> allStrands;

        private static void addStrandsToComplexBuilder(MoleculeDto dto, SimpleConfSpace.Builder builder, List<Strand> strands) {
            if (dto.translateRotate != null) {
                var tr = new StrandFlex.TranslateRotate(
                        dto.translateRotate.rotateDegrees,
                        dto.translateRotate.translateAngstroms
                );

                StreamEx.of(strands).forEach(cs -> builder.addStrand(cs, tr));
            } else {
                StreamEx.of(strands).forEach(builder::addStrand);
            }
        }

        public ForceFieldParamsAndStrands(DesignFileDelegate delegate, AffinityDesign design) {

            forcefieldParams = delegate.frcmodPath == null
                    ? new ForcefieldParams()
                    : new ForcefieldParams(ForcefieldParams.Forcefield.AMBER, new ForcefieldFileParser(
                            getClass().getResourceAsStream(ForcefieldParams.Forcefield.AMBER.paramsPath), Paths.get(delegate.frcmodPath))
            );

            proteinConfSpace = delegate.createConfSpace(design.protein, forcefieldParams);
            ligandConfSpace = delegate.createConfSpace(design.ligand, forcefieldParams);

            allStrands = StreamEx.of(proteinConfSpace.strands, ligandConfSpace.strands)
                    .flatMap(List::stream)
                    .toList();

            var builder = new SimpleConfSpace.Builder();
            addStrandsToComplexBuilder(design.protein, builder, proteinConfSpace.strands);
            addStrandsToComplexBuilder(design.ligand, builder, ligandConfSpace.strands);
            complexConfSpace = builder.build();
        }

        List<Residue> makeAllResidues() {
            return StreamEx.of(allStrands)
                    .flatMap(x -> x.mol.residues.stream())
                    .toList();
        }
    }

    private int makeScanDesigns(AffinityDesign design, List<Residue> allResidues) {
        var target = design.scanSettings.target;
        var residues = design.scanSettings.residues;

        if (target.isEmpty() && residues.isEmpty() || !target.isEmpty() && !residues.isEmpty()) {
            System.err.println("Either target or residues must be specified, but not both");
            return Main.Failure;
        }

        residues.forEach(r -> r.mutability = r.mutability.isEmpty() ? allMutableAaTypes : r.mutability);
        var mutableTargets = residues.isEmpty()
                ? findMutableResiduesAroundTarget(design, scanDistance, target, allResidues)
                : residues;

        return createScanDesigns(design, mutableTargets, scanFlexDistance);
    }

    @NotNull
    private List<ResidueModifier> findMutableResiduesAroundTarget(AffinityDesign design, double dist, String target, List<Residue> allResidues) {
        var targetRes = allResidues.stream()
                .filter(a -> a.getPDBResNumber().equals(target))
                .findFirst()
                .orElseThrow();

        var modifiers = StreamEx.of(allResidues)
                .filter(x -> x != targetRes)
                .filter(x -> !design.scanSettings.excluding.contains(x.getPDBResNumber()))
                .filter(x -> x.distanceTo(targetRes) <= dist)
                .map(CommandBindingAffinity::makeFlexibleResidueModifier)
                .toList();

        modifiers.forEach(a -> a.mutability = allMutableAaTypes);
        return modifiers;
    }

    public static Stream<Residue> nearbyResidues(Residue target, Molecule molecule, double withinDistance) {
        return molecule.residues.stream()
                .filter(x -> x.distanceTo(target) <= withinDistance)
                .filter(x -> target.getChainId() != x.getChainId()
                        || !target.getPDBResNumber().equals(x.getPDBResNumber())
                        || !target.getType().equals(x.getType())
                );
    }

    public static Stream<ResidueModifier> nearbyResidueModifiers(Residue target, Molecule molecule, double withinDistance) {
        return nearbyResidues(target, molecule, withinDistance).map(CommandBindingAffinity::makeFlexibleResidueModifier);
    }

    private int createScanDesigns(AffinityDesign designTemplate, List<ResidueModifier> mutableTargets, double flexDist) {

        var proteinMol = designTemplate.makeProteinMolecule();
        var ligandMol = designTemplate.makeLigandMolecule();

        for (var comboIndices : StreamEx.ofCombinations(mutableTargets.size(), delegate.maxSimultaneousMutations)) {

            var mutableResidues = IntStreamEx.of(comboIndices)
                    .mapToObj(mutableTargets::get)
                    .toList();
            var mutableTargetNums = StreamEx.of(mutableResidues)
                    .map(x -> x.identity.positionIdentifier())
                    .toSet();
            var mutableTargetRes = StreamEx.of(proteinMol.residues)
                    .append(ligandMol.residues)
                    .filter(res -> mutableTargetNums.contains(res.getPDBResNumber()))
                    .toList();

            List<ResidueModifier> proteinResMods = StreamEx.of(mutableTargetRes)
                    .flatMap(residue -> nearbyResidueModifiers(residue, proteinMol, flexDist))
                    .distinct()
                    .toList();

            var ligandResMods = StreamEx.of(mutableTargetRes)
                    .flatMap(residue -> nearbyResidueModifiers(residue, ligandMol, flexDist))
                    .distinct()
                    .toList();

            StreamEx.of(mutableResidues)
                    .zipWith(StreamEx.of(mutableTargetRes))
                    .forKeyValue((residueModifier, residue) -> {
                        if (proteinMol.residues.contains(residue)) {
                            proteinResMods.add(residueModifier);
                        } else {
                            ligandResMods.add(residueModifier);
                        }
                    });

            try {
                var nameTemp = delegate.design.getName().substring(0, delegate.design.getName().lastIndexOf("."));
                var comboName = String.join("-", StreamEx.of(mutableTargetRes)
                        .map(Residue::getPDBResNumber)
                        .toList());

                var path = Path.of(String.format("%s.%s.yaml", nameTemp, comboName));
                Files.deleteIfExists(path);
                var outFile = Files.createFile(path);

                var designCopy = designTemplate.copy();
                designCopy.protein.residueModifiers = StreamEx.of(designCopy.protein.residueModifiers).append(proteinResMods).toList();
                designCopy.ligand.residueModifiers = StreamEx.of(designCopy.ligand.residueModifiers).append(ligandResMods).toList();
                designCopy.scanSettings = null; // no need to copy this into the newly-created design
                designCopy.write(outFile);
            } catch (IOException e) {
                e.printStackTrace();
                return Main.Failure;
            }
        }

        return Main.Success;
    }

    @NotNull
    private static ResidueModifier makeFlexibleResidueModifier(Residue x) {
        var res = new edu.duke.cs.osprey.design.models.Residue();
        res.aminoAcidType = x.getType();
        res.chain = String.valueOf(x.getChainId());
        res.residueNumber = Integer.parseInt(x.getPDBResNumber().substring(1));

        var modifier = new ResidueModifier();
        modifier.identity = res;
        modifier.mutability = List.of();

        return modifier;
    }

    private void printResults(List<ScoredSequence> results) {
        // nop
    }

    private KStarSettings makeKStarSettings(double epsilon) {
        var builder = new KStarSettings.Builder();
        builder.setEpsilon(epsilon);
        builder.addScoreConsoleWriter();
        if (maxNumberConfs > 0) {
            builder.setMaxNumConf(maxNumberConfs);
        }
        builder.setMaxSimultaneousMutations(delegate.maxSimultaneousMutations);
        builder.setStabilityThreshold(stabilityThreshold < 0 ? null : stabilityThreshold);

        if (delegate.writeNConfs > 0) {
            var saveDir = delegate.ensembleDir;
            var scoreWriter = new StructureFileScoreWriter(saveDir, delegate.writeNConfs);
            builder.addScoreWriter(scoreWriter);
        }

        return builder.build();
    }

    private Optional<AffinityDesign> parseAndValidate(File designSpec) {
        return getAffinityDesignFromFile(designSpec);
    }

    @NotNull
    public static Optional<AffinityDesign> getAffinityDesignFromFile(File designSpec) {
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
