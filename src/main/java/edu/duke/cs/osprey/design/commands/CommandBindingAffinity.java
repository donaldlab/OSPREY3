package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.ParametersDelegate;

@Parameters(commandDescription = "Compute a guaranteed-accurate approximation to binding affinity")
public class CommandBindingAffinity implements RunnableCommand {
    public static final String CommandName = "affinity";
    @ParametersDelegate
    private DesignFileDelegate delegate = new DesignFileDelegate();

    @Override
    public int run(JCommander commander, String[] args) {
        return 0;
    }
}
