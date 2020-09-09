package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;

public class StabilityDesign {
    @JsonProperty("osprey_version")
    public String ospreyVersion;

    @JsonProperty("design_name")
    public String designName;

    public float epsilon;

    public MoleculeDto molecule;
}
