package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.gui.io.PDBKt;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;

@Parameters(commandDescription = CommandInvert.CommandDescription)
public class CommandInvert implements CliCommand {

    public static final String CommandName = "invert";
    public static final String CommandDescription = "Flip the sign of the z-coordinate of all atoms in a PDB";

    @Parameter(description = "pdb-file")
    private String pdbFile;

    @Override
    public int run(JCommander commander, String[] args) {
        File f = new File(pdbFile);
        if (!f.exists()) {
            System.err.printf("The file you specified, %s, doesn't exist%n", pdbFile);
            return Main.Failure;
        }


        var molecule = PDBKt.toMolecule(PDBIO.readFile(f), "pdb");
        var invertedCopy = molecule.invertedCopy();
        var invertedStr = PDBKt.toPDB(invertedCopy, false, false, null, null, false, false, false);
        System.out.print(invertedStr);
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
