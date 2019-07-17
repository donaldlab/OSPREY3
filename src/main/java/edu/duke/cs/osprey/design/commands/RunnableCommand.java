package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;

public interface RunnableCommand {
    public int run(JCommander commander, String[] args);
}
