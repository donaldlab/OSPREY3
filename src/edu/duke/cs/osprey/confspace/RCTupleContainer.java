package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Wrapper class for {@link RCTuple} so that we can nicely bundle it with other information,
 * e.g. energy in the case of {@link TupE}.
 */
public class RCTupleContainer {
    public final RCTuple tup;

    RCTupleContainer(RCTuple tup){
        this.tup = tup;
    }

    public static RCTupleContainer fromString(String repr){
        return new RCTupleContainer(RCTuple.fromString(repr));
    }

    @Override
    public String toString(){
        return tup.toString();
    }

}
