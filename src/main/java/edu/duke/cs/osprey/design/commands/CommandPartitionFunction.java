package edu.duke.cs.osprey.design.commands;


import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.ParametersDelegate;
import edu.duke.cs.osprey.design.Main;

@Parameters(commandDescription = CommandPartitionFunction.CommandDescription)
public class CommandPartitionFunction implements RunnableCommand {

    public static final String CommandName = "stability";
    static final String CommandDescription = "Estimate the partition function value(s) of a set of sequence(s)";

    @ParametersDelegate
    private DesignFileDelegate delegate = new DesignFileDelegate();

    private void printHelp(JCommander commander) {
        var msg = String.format("%s: %s", CommandName, CommandDescription);
        System.out.println(msg);
        System.out.println();
        commander.usage();
    }

    private void printMissingDesignFile(JCommander commander) {
        System.out.println("Error: Missing a design file.");
    }

    @Override
    public int run(JCommander commander, String[] args) {
        if (args.length == 1) {
            printHelp(commander);
            return Main.Failure;
        }

        if (delegate.help) {
            printHelp(commander);
            return Main.Success;
        }

        if (delegate.design == null) {
            printMissingDesignFile(commander);
            return Main.Failure;
        }

        return Main.Success;
    }
}
