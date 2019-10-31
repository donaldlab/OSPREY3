package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;

@JsonIgnoreProperties(ignoreUnknown = true)
public class Residue {
    public String chain;
    @JsonProperty("res_num")
    public int residueNumber;
    @JsonProperty("aa_type")
    public AminoAcid aminoAcidType;

    public String positionIdentifier() {
        return String.format("%s%d", chain, residueNumber).toUpperCase();
    }
}
