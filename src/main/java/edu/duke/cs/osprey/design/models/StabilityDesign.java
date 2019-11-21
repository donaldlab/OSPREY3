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

    @JsonProperty("residue_configurations")
    public List<ResidueModifier> residueModifiers;

    public float epsilon;

    public String molecule;

    @JsonProperty("extra_templates")
    public String extraTemplates;

    @JsonProperty("extra_template_coordinates")
    public String extraTemplatesCoordinates;

    public static StabilityDesign parse(File file) throws IOException {
        var mapper = new ObjectMapper(new YAMLFactory());
        var stream = new FileInputStream(file);
        return mapper.readValue(stream, StabilityDesign.class);
    }
}
