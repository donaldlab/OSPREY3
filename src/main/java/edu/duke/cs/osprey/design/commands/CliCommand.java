package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;

public interface CliCommand {
    int run(JCommander commander, String[] args);

    String getCommandName();

    String getCommandDescription();

    default void printHelp(JCommander commander) {
        var msg = String.format("%s: %s", getCommandName(), getCommandDescription());
        System.out.println(msg);
        System.out.println();
        commander.usage();
    }
}
