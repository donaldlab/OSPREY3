package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParametersDelegate;
import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import edu.duke.cs.osprey.design.Main;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Optional;

public abstract class RunnableCommand {

    public abstract int run(JCommander commander, String[] args);

    public abstract String getCommandName();
    public abstract String getCommandDescription();

}
