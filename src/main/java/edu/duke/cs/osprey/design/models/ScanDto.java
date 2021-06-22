package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.List;

public class ScanDto {
    @JsonProperty("target")
    public String target = "";

    @JsonProperty("excluding")
    public List<String> excluding = List.of();

    @JsonProperty("residues")
    public List<ResidueModifier> residues = List.of();
}
