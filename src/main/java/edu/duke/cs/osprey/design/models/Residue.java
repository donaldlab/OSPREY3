package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonProperty;

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
