package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParametersDelegate;
import edu.duke.cs.osprey.design.Main;

import java.util.Optional;

public abstract class RunnableCommand {

    public abstract int run(JCommander commander, String[] args);

    public abstract String getCommandName();
    public abstract String getCommandDescription();

    @ParametersDelegate
    private DesignFileDelegate delegate = new DesignFileDelegate();

    Optional<Integer> processHelpAndNoArgs(JCommander commander, String[] args) {
        if (args.length == 1) {
            printHelp(commander);
            return Optional.of(Main.Failure);
        }

        if (delegate.help) {
            printHelp(commander);
            return Optional.of(Main.Success);
        }

        if (delegate.design == null) {
            printMissingDesignFile(commander);
            return Optional.of(Main.Failure);
        }

        return Optional.empty();
    }

    void printHelp(JCommander commander) {
        var msg = String.format("%s: %s", getCommandName(), getCommandDescription());
        System.out.println(msg);
        System.out.println();
        commander.usage();
    }

    void printMissingDesignFile(JCommander commander) {
        System.out.println("Error: Missing a design file.");
    }
}
