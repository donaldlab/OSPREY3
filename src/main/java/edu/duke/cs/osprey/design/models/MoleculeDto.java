package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.List;

public class MoleculeDto {
    @JsonProperty("residue_configurations")
    public List<ResidueModifier> residueModifiers = List.of();

    @JsonProperty("extra_templates")
    public String extraTemplates = "";

    @JsonProperty("extra_templates_coordinates")
    public String extraTemplatesCoordinates = "";

    @JsonProperty("coordinates")
    public String coordinates;

    @JsonProperty("extra_rotamers")
    public String additionalRotamers = "";
}
