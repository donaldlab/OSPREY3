package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.internal.Lists;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Sets;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.models.AffinityDesign;
import edu.duke.cs.osprey.design.models.ResidueModifier;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Parameters(commandDescription = CommandMakeFlexShell.CommandDescription)
public class CommandMakeFlexShell extends RunDesignCommand {
    public static final String CommandName = "add-flexible-shell";
    public static final String CommandDescription = "Adds a flexible shell around mutable residues";

    @Parameter(names = "--flex-distance", description = "Distance (in angstroms) around mutable residues in which residues are flexible.")
    public double flexDistance = 3.0;

    @Override
    public int run(JCommander commander, String[] args) {
        return processHelpAndNoArgs(commander, args)
                .orElseGet(() -> parseAndValidate(delegate.design)
                        .map(this::addFlexibleShell)
                        .orElse(Main.Failure));
    }

    private Stream<ResidueModifier> mutableResidues(AffinityDesign design) {
        return Stream.of(design.ligand.residueModifiers, design.protein.residueModifiers)
                .flatMap(Collection::stream)
                .filter(x -> !x.mutability.isEmpty());
    }

    private List<ResidueModifier> makeResidueModifiersInMolecule(AffinityDesign design, Molecule molecule) {
        var mutableResidues = mutableResidues(design).collect(Collectors.toSet());

        return mutableResidues(design)
                .map(mod -> toMoleculeResidue(List.of(design.makeProteinMolecule(), design.makeLigandMolecule()), mod))
                .flatMap(residue -> CommandBindingAffinity.nearbyResidueModifiers(residue, molecule, flexDistance))
                .reduce(ImmutableList.of(),
                        (intermediateList, residueModifier) -> {
                            var residues = Stream.concat(intermediateList.stream().map(mod -> mod.identity),
                                    mutableResidues(design).map(mod -> mod.identity)).collect(Collectors.toList());

                            return residues.contains(residueModifier.identity)
                                    ? intermediateList
                                    : new ImmutableList.Builder<ResidueModifier>().addAll(intermediateList).add(residueModifier).build();
                        },
                        (intermediate, toAppend) -> ImmutableList.copyOf(Iterables.concat(intermediate, toAppend)));
    }

    private int addFlexibleShell(AffinityDesign design) {
        var protein = design.makeProteinMolecule();
        var ligand = design.makeLigandMolecule();

        var proteinModifiers = makeResidueModifiersInMolecule(design, protein);
        var ligandModifiers = makeResidueModifiersInMolecule(design, ligand);
        design.protein.residueModifiers.addAll(proteinModifiers);
        design.ligand.residueModifiers.addAll(ligandModifiers);

        try {
            var nameTemp = delegate.design.getName().substring(0, delegate.design.getName().lastIndexOf("."));
            var path = Path.of(String.format("%s.flex-shell.yaml", nameTemp)).getFileName();
            design.write(path);
        } catch (IOException e) {
            e.printStackTrace();
        }

        return Main.Success;
    }

    private boolean isSamePosition(Residue residue, ResidueModifier modifier) {
        return residue.getChainId() == modifier.identity.chain.charAt(0)
                && modifier.identity.aminoAcidType.equals(residue.getType())
                && modifier.identity.residueNumber == Integer.parseInt(residue.getPDBResNumber().substring(1));
    }

    private Residue toMoleculeResidue(Collection<Molecule> molecule, ResidueModifier mod) {
        return molecule.stream()
                .flatMap(mol -> mol.residues.stream())
                .filter(residue -> isSamePosition(residue, mod))
                .findAny()
                .orElseThrow();
    }

    private Optional<AffinityDesign> checkFlexShellRequirements(AffinityDesign design) {
        var hasMutableResidue = mutableResidues(design).anyMatch(mod -> !mod.mutability.isEmpty());

        if (!hasMutableResidue) {
            System.err.print("Making a flex shell requires at least one mutable residue to flex around%n");
            return Optional.empty();
        }

        return Optional.of(design);
    }

    private Optional<AffinityDesign> parseAndValidate(File designFile) {
        var design = CommandBindingAffinity.getAffinityDesignFromFile(designFile);
        return design.flatMap(this::checkFlexShellRequirements);
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
