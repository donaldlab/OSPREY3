package edu.duke.cs.osprey.design.commands;


import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.ParametersDelegate;
import edu.duke.cs.osprey.design.Main;
import org.jetbrains.annotations.Nullable;

import java.util.Optional;

@Parameters(commandDescription = CommandPartitionFunction.CommandDescription)
public class CommandPartitionFunction extends RunnableCommand {

    public static final String CommandName = "stability";
    static final String CommandDescription = "Estimate the partition function value(s) of a set of sequence(s)";

    @Override
    public int run(JCommander commander, String[] args) {
        var retVal = processHelpAndNoArgs(commander, args);

        if (retVal.isPresent()) {
            return retVal.get();
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
