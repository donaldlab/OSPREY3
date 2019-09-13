package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.converters.FileConverter;

import java.io.File;

public class DesignFileDelegate {

    @Parameter(description = "Path to design file. For more information on the syntax of the design file, use --help-design", names = {"--design", "-d"}, converter = FileConverter.class, validateWith=FileExistsValidation.class)
    public File design;

    @Parameter(description = "Prints this help information", names={"--help", "-h"}, help = true)
    boolean help;

    @Parameter(description = "Prints this help information", names="--help-design")
    boolean helpDesign;

    @Parameter(description = "The approximation accuracy. Z* = (1 - epsilon)Z", names={"--epsilon", "-e"})
    double epsilon;
}
