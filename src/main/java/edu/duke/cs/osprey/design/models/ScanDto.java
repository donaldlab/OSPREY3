package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.List;

public class ScanDto {
    @JsonProperty("target")
    public String target = "";

    @JsonProperty("distance")
    public double distance;

    @JsonProperty("excluding")
    public List<String> excluding = List.of();
}
