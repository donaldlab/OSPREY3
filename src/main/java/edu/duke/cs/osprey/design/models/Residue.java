package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.Objects;

@JsonIgnoreProperties(ignoreUnknown = true)
public class Residue {

    public String chain;

    @JsonProperty("res_num")
    public int residueNumber;

    @JsonProperty("aa_type")
    public String aminoAcidType; // 3-letter-code

    public String positionIdentifier() {
        return String.format("%s%d", chain, residueNumber).toUpperCase();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Residue residue = (Residue) o;
        return residueNumber == residue.residueNumber &&
                chain.equals(residue.chain) &&
                aminoAcidType.equals(residue.aminoAcidType);
    }

    @Override
    public int hashCode() {
        return Objects.hash(chain, residueNumber, aminoAcidType);
    }
}
