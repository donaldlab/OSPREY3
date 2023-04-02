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

    @ParametersDelegate
    protected DesignFileDelegate delegate = new DesignFileDelegate();

    Optional<Integer> processHelpAndNoArgs(JCommander commander, String[] args) {
        if (args.length == 1) {
            printHelp(commander);
            return Optional.of(Main.Failure);
        }

        if (delegate.help) {
            printHelp(commander);
            return Optional.of(Main.Success);
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

    public <T> Optional<T> parseDesignSpec(Class<T> t) {
        try {
            var mapper = new ObjectMapper(new YAMLFactory());
            mapper.configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, false);
            var stream = new FileInputStream(delegate.design);
            return Optional.of(mapper.readValue(stream, t));
        } catch (IOException e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }
}
