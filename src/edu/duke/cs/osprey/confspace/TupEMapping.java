package edu.duke.cs.osprey.confspace;

public class TupEMapping extends TupE {
    public RCTuple mappedTup;

    public TupEMapping(RCTuple tup, double E) {
        super(tup, E);
        this.mappedTup = null;
    }

    public TupEMapping(RCTuple tup, RCTuple mapTup, double E){
        super(tup, E);
        this.mappedTup = mapTup;
    }

    public static TupEMapping fromString(String repr) {
        String[] tupSplits = repr.split("->");
        TupE tupE = TupE.fromString(tupSplits[0]+"->"+tupSplits[2]);

        return new TupEMapping(tupE.tup, RCTuple.fromString(tupSplits[1]), tupE.E);
    }

    @Override
    public String toString() {
        return tup.toString()+"->"+mappedTup.stringListing()+"->"+E;
    }

    public String toString_long() {
        return tup.stringListing()+"->"+mappedTup.stringListing()+"->"+E;
    }


}
