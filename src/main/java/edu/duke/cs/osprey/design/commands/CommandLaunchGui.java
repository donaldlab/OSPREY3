package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.gui.MainKt;

import edu.duke.cs.osprey.design.Main;

@Parameters(commandDescription = CommandLaunchGui.CommandDescription)
public class CommandLaunchGui implements CliCommand {

    public static final String CommandName = "setup-design";
    public static final String CommandDescription = "Prep structures and specify conformation spaces in a graphical user application";

    @Override
    public int run(JCommander commander, String[] args) {
        MainKt.main();
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
