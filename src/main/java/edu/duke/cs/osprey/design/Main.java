package edu.duke.cs.osprey.design;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.MissingCommandException;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import edu.duke.cs.osprey.design.commands.*;

import java.util.Map;

import static java.lang.System.exit;

public class Main {

    @Parameter(names = {"--help", "-h"}, help = true)
    private boolean help;

    public final static int Success = 0;
    public final static int Failure = 1;
    public final static String ProgramName = "osprey";

    private void printErrorMessage(String msg) {
        System.out.println(msg);
        System.out.println();
        System.out.println("Use --help or <command> --help for more info.");
        System.out.println();
    }

    public int run(String[] args) {

        var commandMap = Map.of(
                CommandLaunchGui.CommandName, new CommandLaunchGui(),
//                CommandPartitionFunction.CommandName, new CommandPartitionFunction(),
//                CommandBindingAffinity.CommandName, new CommandBindingAffinity(),
//                CommandTopNConfs.CommandName, new CommandTopNConfs(),
//                CommandGMEC.CommandName, new CommandGMEC(),
//                CommandMakeFlexShell.CommandName, new CommandMakeFlexShell(),
//                CommandKStar.CommandName, new CommandKStar(),
                CompiledConfSpaceKStar.CommandName, new CompiledConfSpaceKStar(),
                CompiledConfSpaceBBKStar.CommandName, new CompiledConfSpaceBBKStar(),
                CommandInvert.CommandName, new CommandInvert()
        );

        var builder = JCommander.newBuilder()
                .programName(ProgramName)
                .addObject(this);
        commandMap.forEach(builder::addCommand);
        var commander = builder.build();

        try {
            commander.parse(args);
        } catch (MissingCommandException ex) {
            printErrorMessage(String.format("Error: the command \"%s\" does not exist.", ex.getUnknownCommand()));
            commander.usage();
            return Failure;
        } catch (ParameterException ex) {
            printErrorMessage(String.format("Error: %s", ex.getMessage()));
            return Failure;
        }

        if (args.length == 0 || this.help) {
            commander.usage();
            return Success;
        }

        var parsedCommand = commander.getParsedCommand();
        return commandMap.get(parsedCommand).run(commander.getCommands().get(parsedCommand), args);
    }

    public static void main(String[] args) {

        Thread.setDefaultUncaughtExceptionHandler((thread, throwable) -> {
            if (throwable instanceof OutOfMemoryError) {
                printMemoryStatistics();
            }

            throwable.printStackTrace();
            exit(Failure);
        });

        exit(new Main().run(args));
    }

    private static void printMemoryStatistics() {
        var runtime = Runtime.getRuntime();
        System.err.printf(
                "Current HeapSize: %d\nMax Heap Size: %d\nCurrent Free Space: %d\n%n",
                runtime.totalMemory(),
                runtime.maxMemory(),
                runtime.freeMemory()
        );
    }
}

