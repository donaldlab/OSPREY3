package edu.duke.cs.osprey.kstar;

public class InitException extends RuntimeException {

    public InitException(ConfSpaceType type, String name) {
        super(String.format("set %s for the %s conf space info before running", name, type.name()));
    }
}
