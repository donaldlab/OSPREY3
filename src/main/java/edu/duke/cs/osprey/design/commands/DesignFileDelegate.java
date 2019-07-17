package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.converters.FileConverter;

import java.io.File;

public class DesignFileDelegate {

    @Parameter(description = "Path to design file", names = "--design", converter = FileConverter.class, validateWith=FileExistsValidation.class)
    public File design;

    @Parameter(names={"--help", "-h"}, help = true)
    boolean help;
}
