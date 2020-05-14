package edu.duke.cs.osprey.confspace;

public class TupMapping extends RCTupleContainer{
    public RCTuple mapTup;

    public TupMapping(RCTuple tup) {
        super(tup);
    }

    public TupMapping(RCTuple tup, RCTuple mapTup) {
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
