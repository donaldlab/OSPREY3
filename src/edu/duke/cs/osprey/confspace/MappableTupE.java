package edu.duke.cs.osprey.confspace;

import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MappableTupE extends TupE {
    public RCTuple mappedTup;

    public MappableTupE(RCTuple tup, double E) {
        super(tup, E);
        this.mappedTup = null;
    }

    public MappableTupE(RCTuple tup, RCTuple mapTup, double E){
        super(tup, E);
        this.mappedTup = mapTup;
    }

    public MappableTupE(@NotNull String repr) {
        super(repr.split("->")[0]+"->"+repr.split("->")[2]);

        ArrayList<Integer> mapPos = new ArrayList<>();
        ArrayList <Integer> mapRCs = new ArrayList<>();

        // Form arrays from mapTupString
        Pattern point = Pattern.compile("\\d+=\\d+");
        Matcher m = point.matcher(repr.split("->")[1]);
        while (m.find()){
            String[] splits = m.group().split("=");
            mapPos.add(Integer.parseInt(splits[0]));
            mapRCs.add(Integer.parseInt(splits[1]));
        }
        this.mappedTup = new RCTuple(mapPos, mapRCs);
    }

    @Override
    public String toString() {
        return tup.stringListing()+"->"+mappedTup.stringListing()+"->"+E;
    }

    @Override
    public String toString_short() {
        return tup.toString()+"->"+mappedTup.toString()+"->"+E;
    }


}
