package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;
import edu.duke.cs.osprey.design.*;
import edu.duke.cs.osprey.gui.OspreyGui;
import edu.duke.cs.osprey.gui.io.ConfLib;
import edu.duke.cs.osprey.gui.io.ConfLibKt;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.PDBIO;
import org.joml.Vector3d;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

@Parameters(commandDescription = CommandFlipChirality.CommandDescription)
public class CommandFlipChirality extends RunnableCommand {

    public static final String CommandName = "flip-chirality";
    public static final String CommandDescription = "Reflects a molecule about its center of mass in the Z plane.";

    @Parameter(names = "--pdb", description = "The path to the PDB file to reflect.",
            converter = FileConverter.class,
            validateWith = FileExistsValidation.class)
    public File pdb;

    private String[] args;

    private double[] calcCenterOfMass(List<Atom> atoms) {
        var xMass = atoms.stream().map(atom -> atom.getCoords()[0]).mapToDouble(f -> f).sum();
        var yMass = atoms.stream().map(atom -> atom.getCoords()[1]).mapToDouble(f -> f).sum();
        var zMass = atoms.stream().map(atom -> atom.getCoords()[2]).mapToDouble(f -> f).sum();

        return new double[] { xMass / atoms.size(), yMass / atoms.size(), zMass / atoms.size()};
    }

    @Override
    public int run(JCommander commander, String[] args) {
        this.args = args;
        var molecule = PDBIO.readFile(pdb);
        var normal = new double [] {0, 0, 1};

        var centerOfMass =
                calcCenterOfMass(molecule.residues.stream().flatMap(x -> x.atoms.stream()).toList());

        for (var residue : molecule.residues) {
            for (var atom : residue.atoms) {
                var coords = atom.getCoords();
                atom.setCoords(coords[0], coords[1], -coords[2]);
            }
        }

        var confLibStr = OspreyGui.INSTANCE.getResourceAsString("conflib/lovell.conflib", Charset.defaultCharset());
        var confLib = ConfLib.Companion.from(confLibStr);

        for (var fragment: confLib.getFragments().entrySet()) {
            var frag = fragment.getValue();
            for (var conf: frag.getConfs().entrySet()) {
                var c = conf.getValue();
                for (Map.Entry<ConfLib.AtomInfo, Vector3d> atomInfoVector3dEntry : c.getCoords().entrySet()) {
                    atomInfoVector3dEntry.getValue().z = - atomInfoVector3dEntry.getValue().z;
                }

                for (Map.Entry<ConfLib.Anchor, ConfLib.AnchorCoords> anchorAnchorCoordsEntry : c.getAnchorCoords().entrySet()) {
                    var ac = anchorAnchorCoordsEntry.getValue();
                    ac.getCoords(0).z = -ac.getCoords(0).z;
                    ac.getCoords(1).z = -ac.getCoords(1).z;
                    ac.getCoords(2).z = -ac.getCoords(2).z;
                }
            }
        }
        var toml = ConfLibKt.toToml(confLib, null);
        var newPath = pdb.toPath() + ".D.pdb";
        var tomlPath = pdb.toPath() + ".toml";
        try {
            Files.writeString(Path.of(tomlPath), toml);
        } catch (IOException e) {
            e.printStackTrace();
        }
        PDBIO.writeFile(molecule, newPath);
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
