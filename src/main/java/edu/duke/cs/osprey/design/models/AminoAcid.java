package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonIdentityInfo;
import com.fasterxml.jackson.annotation.JsonValue;

public enum AminoAcid {

    Alanine("Ala", 'A'),
    Arginine("Arg", 'R'),
    Asparagine("Asn", 'N'),
    Aspartate("Asp", 'D'),
    Cysteine("Cys", 'C'),
    Glutamine("Gln", 'Q'),
    Glutamate("Glu", 'E'),
    Glycine("Gly", 'G'),
    Histidine("His", 'H'),
    HistidineD("Hid", 'D'),
    HistidineE("Hie", 'O'),
    HistidineP("Hip", 'Z'),
    Isoleucine("Ile", 'I'),
    Leucine("Leu", 'L'),
    Lysine("Lys", 'K'),
    Methionine("Met", 'M'),
    Phenylalanine("Phe", 'F'),
    Proline("Pro", 'P'),
    Serine("Ser", 'S'),
    Threonine("Thr", 'T'),
    Tryptophan("Trp", 'W'),
    Tyrosine("Tyr", 'Y'),
    Valine("Val", 'V'),
    Termination("TERM", (char)0);


    private final String threeLetter;
    private final char oneLetter;

    AminoAcid(String threeLetter, char oneLetter) {
        this.threeLetter = threeLetter;
        this.oneLetter = oneLetter;
    }

    @JsonValue
    public String toValue() {
        return threeLetter;
    }

}


