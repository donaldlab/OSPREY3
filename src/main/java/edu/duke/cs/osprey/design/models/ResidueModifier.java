package edu.duke.cs.osprey.design.models;

import java.util.List;
import java.util.Objects;
import java.util.Set;

public class ResidueModifier {
    public List<String> mutability = List.of();
    public Flexibility flexibility = new Flexibility(); // default flexibility
    public Residue identity;


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ResidueModifier that = (ResidueModifier) o;
        return mutability.equals(that.mutability) &&
                flexibility.equals(that.flexibility) &&
                identity.equals(that.identity);
    }

    @Override
    public int hashCode() {
        return Objects.hash(mutability, flexibility, identity);
    }
}
