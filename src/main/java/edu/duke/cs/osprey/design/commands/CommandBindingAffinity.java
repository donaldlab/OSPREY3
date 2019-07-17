package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.ParametersDelegate;
import edu.duke.cs.osprey.design.Main;

@Parameters(commandDescription = CommandBindingAffinity.CommandDescription)
public class CommandBindingAffinity extends RunnableCommand {

    public static final String CommandName = "affinity";
    public static final String CommandDescription = "Compute a guaranteed-accurate approximation to binding affinity (K*).";

    @Override
    public int run(JCommander commander, String[] args) {
        var opt = processHelpAndNoArgs(commander, args);
        if (opt.isPresent()) {
            return opt.get();
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
