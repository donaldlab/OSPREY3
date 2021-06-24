package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.IParameterValidator;
import com.beust.jcommander.ParameterException;

import java.io.File;

public class FileExistsValidation implements IParameterValidator {

    @Override
    public void validate(String name, String value) throws ParameterException {
        var file = new File(value);
        if (!file.exists()) {
            throw new ParameterException("File " + value + " does not exist");
        }
    }
}
