package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.models.MoleculeDesign;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.structure.PDBIO;

import java.nio.file.Paths;

@Parameters(commandDescription = CommandGMEC.CommandDescription)
public class CommandGMEC extends RunnableCommand {

    public static final String CommandName = "gmec";
    static final String CommandDescription = "Find the GMEC of the conformation space";

    @Override
    public int run(JCommander commander, String[] args) {

        var retVal = processHelpAndNoArgs(commander, args);

        if (retVal.isPresent()) {
            return retVal.get();
        }

        var designOpt = parseDesignSpec(MoleculeDesign.class);
        if (designOpt.isEmpty()) {
            return Main.Failure;
        }

        return runGMEC(designOpt.get());
    }

    @Override
    public String getCommandName() {
        return CommandName;
    }

    @Override
    public String getCommandDescription() {
        return CommandDescription;
    }

    private int runGMEC(MoleculeDesign design) {
        var ffParams = new ForcefieldParams();
        var confSpace = delegate.createConfSpace(design.molecule, ffParams);

        var eCalc = new EnergyCalculator.Builder(confSpace, ffParams)
                .setParallelism(delegate.getParallelism())
                .build();

        var confEnergyCalc = new ConfEnergyCalculator.Builder(confSpace, eCalc)
                .build();

        var emat = new SimplerEnergyMatrixCalculator.Builder(confEnergyCalc)
                .build()
                .calcEnergyMatrix();

        var astar = new ConfAStarTree.Builder(emat, confSpace)
                .setMPLP()
                .build();

        var gmec = new SimpleGMECFinder.Builder(astar, confEnergyCalc)
                .build()
                .find();

        var structure = confSpace.makeMolecule(gmec.getAssignments());

        if (delegate.saveDir != null) {
            var outFile = Paths.get(delegate.saveDir, "gmec.pdb").toFile();
            PDBIO.writeFile(structure.mol, "GMEC structure", gmec.getEnergy(), outFile);
        }

        return Main.Success;
    }
}
