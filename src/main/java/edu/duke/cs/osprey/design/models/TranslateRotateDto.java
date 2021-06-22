package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import edu.duke.cs.osprey.confspace.StrandFlex;

public class TranslateRotateDto {
    @JsonProperty("degrees_rotation")
    public double rotateDegrees = StrandFlex.TranslateRotate.DefaultMaxRotDegrees;

    @JsonProperty("angstroms_translation")
    public double translateAngstroms = StrandFlex.TranslateRotate.DefaultMaxTranslation;

    @JsonCreator
    public TranslateRotateDto(@JsonProperty("degrees_rotation") double rotateDegrees, @JsonProperty("angstroms_translation") double translateAngstroms) {
        this.rotateDegrees = rotateDegrees;
        this.translateAngstroms = translateAngstroms;
    }

    @JsonCreator
    public TranslateRotateDto(boolean use) {
    }
}
