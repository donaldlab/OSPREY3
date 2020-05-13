package edu.duke.cs.osprey.confspace;

public class TupMapping extends RCTupleContainer{
    public RCTuple mapTup;

    TupMapping(RCTuple tup) {
        super(tup);
    }

    TupMapping(RCTuple tup, RCTuple mapTup) {
        super(tup);
        this.mapTup = mapTup;
    }

    @Override
    public String toString(){
        return String.join("->",tup.toString(), mapTup.toString());
    }

    public static TupMapping fromString(String repr){
        String[] parts = repr.split("->");
        return new TupMapping(RCTuple.fromString(parts[0]), RCTuple.fromString(parts[0]));
    }
}
