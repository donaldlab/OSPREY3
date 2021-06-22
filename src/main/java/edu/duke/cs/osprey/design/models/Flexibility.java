package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.Objects;

public class Flexibility {
    @JsonProperty("is_flexible")
    public boolean isFlexible = true;
    @JsonProperty("include_structure_rotamer")
    public boolean includeStructureRotamer = true;
    @JsonProperty("use_continuous")
    public boolean useContinuous = true;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Flexibility that = (Flexibility) o;
        return isFlexible == that.isFlexible &&
                includeStructureRotamer == that.includeStructureRotamer &&
                useContinuous == that.useContinuous;
    }

    @Override
    public int hashCode() {
        return Objects.hash(isFlexible, includeStructureRotamer, useContinuous);
    }
}
