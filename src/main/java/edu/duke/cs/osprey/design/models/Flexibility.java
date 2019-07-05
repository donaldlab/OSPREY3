package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonProperty;

public class Flexibility {
    @JsonProperty("is_flexible")
    public boolean isFlexible;
    @JsonProperty("include_structure_rotamer")
    public boolean includeStructureRotamer;
}
