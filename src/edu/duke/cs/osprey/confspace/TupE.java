/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * Like ConfPair but with an (RCTuple,energy) pair
 * 
 * @author mhall44
 */
public class TupE implements Comparable<TupE> {
    
    public RCTuple tup;
    public double E;

    public TupE(RCTuple tup, double E) {
        this.tup = tup;
        this.E = E;
    }
    public TupE(String repr){
        String[] parts = repr.split("->");
        String tupString = parts[0];
        String eString = parts[1];

        ArrayList <Integer> pos = new ArrayList<Integer>();
        ArrayList <Integer> RCs = new ArrayList<Integer>();

        // Form arrays from tupString
        Pattern point = Pattern.compile("\\d+=\\d+");
        Matcher m = point.matcher(tupString);
        while (m.find()){
            String[] splits = m.group().split("=");
            pos.add(Integer.parseInt(splits[0]));
            RCs.add(Integer.parseInt(splits[1]));
        }

        this.tup = new RCTuple(pos, RCs);
        this.E = Double.parseDouble(eString);
    }

    public TupE permute(int[] perm){
        return new TupE(this.tup.permutedCopy(perm), this.E);
    }
    
	
    @Override
    public int compareTo(TupE o) {
        return Double.compare(E,o.E);
    }

    @Override
    public String toString() {
        return tup.stringListing()+"->"+E;
    }

    public String toString_short() {
        return tup.toString()+"->"+E;
    }

    public boolean equals(TupE o ){
        boolean value = false;
        if (E == o.E){
            if (tup.equals(o.tup)){
                value =  true;
            }
        }
        return value;
    }

    public int size(){
        return this.tup.size();
    }
    
}
