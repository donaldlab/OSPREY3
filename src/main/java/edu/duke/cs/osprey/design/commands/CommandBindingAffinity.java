package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.models.AffinityDesign;
import edu.duke.cs.osprey.design.models.Flexibility;
import edu.duke.cs.osprey.design.models.ResidueModifier;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import org.eclipse.collections.impl.factory.Sets;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

@Parameters(commandDescription = CommandBindingAffinity.CommandDescription)
public class CommandBindingAffinity extends RunnableCommand {

    public static final String CommandName = "affinity";
    public static final String CommandDescription = "Compute an epsilon approximation to binding affinity (K*).";

    @Parameter(names = "--do-scan", description = "Runs a scan using the scan settings specified in the design.")
    public boolean doScan;

    @Parameter(names = "--scan-output", description = "Specifies the output directory to save the scan designs in.")
    public String scanOutput;

    @Override
    public int run(JCommander commander, String[] args) {
        var opt = processHelpAndNoArgs(commander, args);
        if (opt.isPresent()) {
            return opt.get();
        }

        return parseAndValidate(delegate.design)
                .map(this::runAffinityDesign)
                .orElse(Main.Failure);
    }

    private int runAffinityDesign(AffinityDesign design) {
        var forcefieldParams = new ForcefieldParams();
        var confSpace1 = delegate.createConfSpace(design.protein, forcefieldParams);
        var confSpace2 = delegate.createConfSpace(design.ligand, forcefieldParams);
        var strands = List.of(confSpace1.strands, confSpace2.strands)
                .stream()
                .flatMap(List::stream)
                .collect(Collectors.toList());

        var complexConfSpace = new SimpleConfSpace.Builder()
                .addStrands(strands)
                .build();

        if (doScan && design.scanSettings != null) {
            var dist = design.scanSettings.distance;
            var target = design.scanSettings.target;
            var residues = design.scanSettings.residues;

            if (target.isEmpty() && residues.isEmpty() || !target.isEmpty() && !residues.isEmpty()) {
                System.err.println("Either target or residues must be specified, but not both");
                return Main.Failure;
            }

            var allResidues = strands.stream().flatMap(x -> x.mol.residues.stream()).collect(Collectors.toList()) ;

            var mutableTargets = residues.isEmpty()
                    ? findMutableResiduesAroundTarget(design, dist, target, allResidues)
                    : specifyMutableResiduesInDesign(allResidues, design.scanSettings.excluding, residues);

            return createScanDesigns(design, mutableTargets, allResidues, 4);
        }

        // Exit early if just trying to validate input
        if (delegate.verifyInput) {
            System.out.println("Design file validated.");
            return Main.Success;
        }

        /* Decides whether to use CPU(s) and/or GPU(s) (purely implementation specific) */
        var parallelism = delegate.getParallelism();

        /* Used to calculate energies of a molecule, also used to minimize the molecule */
        var minimizingECalc = new EnergyCalculator.Builder(complexConfSpace, forcefieldParams)
                .setParallelism(parallelism)
                .build();

        var epsilon = delegate.epsilon > 0 ? delegate.epsilon : design.epsilon;
        var kstar = new KStar(confSpace1, confSpace2, complexConfSpace, makeKStarSettings(epsilon));

        for (var info : kstar.confSpaceInfos()) {
            var referenceEnergies = new SimpleReferenceEnergies.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc).build();

            info.confEcalc = new ConfEnergyCalculator.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc)
                    .setReferenceEnergies(referenceEnergies)
                    .build();

            var minimizedEnergyMatrix = new SimplerEnergyMatrixCalculator.Builder(((SimpleConfSpace) info.confSpace), minimizingECalc)
                    .build()
                    .calcEnergyMatrix();

            info.pfuncFactory = (rcs) -> new GradientDescentPfunc(
                    info.confEcalc,
                    new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                    new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                    rcs.getNumConformations()
            );
        }

        printResults(kstar.run(minimizingECalc.tasks));
        return Main.Success;
    }

    private List<Residue> specifyMutableResiduesInDesign(List<Residue> allResidues, List<String> excludingResidues, List<String> specifiedResidues) {
        return allResidues.stream()
                .filter(x -> !excludingResidues.contains(x.getPDBResNumber()))
                .filter(x -> specifiedResidues.contains(x.getPDBResNumber()))
                .collect(Collectors.toList());
    }

    @NotNull
    private List<Residue> findMutableResiduesAroundTarget(AffinityDesign design, double dist, String target, List<Residue> allResidues) {
        var targetRes = allResidues.stream()
                .filter(a -> a.getPDBResNumber().equals(target))
                .findFirst()
                .orElseThrow();

        return allResidues.stream()
                .filter(x -> x != targetRes)
                .filter(x -> !design.scanSettings.excluding.contains(x.getPDBResNumber()))
                .filter(x -> x.distanceTo(targetRes) <= dist)
                .collect(Collectors.toList());
    }

    private int createScanDesigns(AffinityDesign designTemplate, List<Residue> mutableTargets, List<Residue> allResidues, double flexDist) {
        for (var residue : mutableTargets) {

            var flexShellResidues = allResidues.stream()
                    .filter(x -> x.distanceTo(residue) <= flexDist)
                    .filter(x -> !residue.equals(x))
                    .map(Residue::getPDBResNumber)
                    .collect(Collectors.toSet());

            var proteinMol = PDBIO.read(designTemplate.protein.coordinates);
            var ligandMol = PDBIO.read(designTemplate.ligand.coordinates);

            var proteinRes = proteinMol.residues.stream()
                    .map(Residue::getPDBResNumber)
                    .collect(Collectors.toSet());

            var ligandRes = ligandMol.residues.stream()
                    .map(Residue::getPDBResNumber)
                    .collect(Collectors.toSet());

            var proteinFlexStr = Sets.intersect(flexShellResidues, proteinRes);
            var ligandFlexStr = Sets.intersect(flexShellResidues, ligandRes);

            var proteinFlex = proteinMol.residues.stream()
                    .filter(x -> proteinFlexStr.contains(x.getPDBResNumber()))
                    .collect(Collectors.toUnmodifiableList());

            var ligandFlex = ligandMol.residues.stream()
                    .filter(x -> ligandFlexStr.contains(x.getPDBResNumber()))
                    .collect(Collectors.toUnmodifiableList());

            var proteinResMods = proteinFlex.stream()
                    .map(this::makeFlexibleResidueModifier)
                    .collect(Collectors.toList());

            var ligandResMods = ligandFlex.stream()
                    .map(this::makeFlexibleResidueModifier)
                    .collect(Collectors.toList());

            var mutableResModifier = makeFlexibleResidueModifier(residue);
            mutableResModifier.mutability = List.of("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
                    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL");

            if (proteinMol.residues.contains(residue)) {
                proteinResMods.add(mutableResModifier);
            } else {
                ligandResMods.add(mutableResModifier);
            }

            try {
                var nameTemp = delegate.design.getName().substring(0, delegate.design.getName().lastIndexOf("."));
                var path = Path.of(String.format("%s.%s.yaml", nameTemp, residue.getPDBResNumber()));
                Files.deleteIfExists(path);
                var outFile = Files.createFile(path);

                var designCopy = designTemplate.copy();
                designCopy.protein.residueModifiers = proteinResMods;
                designCopy.ligand.residueModifiers = ligandResMods;
                designCopy.write(outFile);
            } catch (IOException e) {
                e.printStackTrace();
                return Main.Failure;
            }
        }

        return Main.Success;
    }

    @NotNull
    private ResidueModifier makeFlexibleResidueModifier(Residue x) {
        var res = new edu.duke.cs.osprey.design.models.Residue();
        res.aminoAcidType = x.getType();
        res.chain = String.valueOf(x.getChainId());
        res.residueNumber = Integer.parseInt(x.getPDBResNumber().substring(1));

        var flexibility = new Flexibility();
        flexibility.includeStructureRotamer = true;
        flexibility.isFlexible = true;
        flexibility.useContinuous = true;

        var modifier = new ResidueModifier();
        modifier.identity = res;
        modifier.flexibility = flexibility;
        modifier.mutability = List.of();

        return modifier;
    }

    private void printResults(List<KStar.ScoredSequence> results) {
        for (var result : results) {
            System.out.println("result: ");
            System.out.printf("\tsequence: %s%n", result.sequence);
            System.out.printf("\tscore: %s%n", result.score);
        }
    }

    private KStar.Settings makeKStarSettings(double epsilon) {
        return new KStar.Settings.Builder()
                .setEpsilon(epsilon)
                .addScoreConsoleWriter()
                .build();
    }

    private BBKStar.Settings makeBBKStarSettings() {
        return new BBKStar.Settings.Builder().build();
    }

    static Optional<AffinityDesign> parseAndValidate(File designSpec) {

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
