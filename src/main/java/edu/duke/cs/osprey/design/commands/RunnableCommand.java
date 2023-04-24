package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParametersDelegate;
import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import edu.duke.cs.osprey.design.Main;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import static jogamp.common.os.elf.SectionArmAttributes.Tag.File;

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

    void cleanupStuffFromPreviousRuns() {
        // Reading only files in the directory
        try {
            List<java.io.File> filesInPwd = Files.list(Paths.get(""))
                    .map(Path::toFile)
                    .filter(java.io.File::isFile)
                    .filter(f -> f.toString().contains(".confdb") || f.toString().contains(".wal"))
                    .toList();

            filesInPwd.forEach(java.io.File::delete);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void printHelp(JCommander commander) {
        var msg = String.format("%s: %s", getCommandName(), getCommandDescription());
        System.out.println(msg);
        System.out.println();
        commander.usage();
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
