/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.structure;

import java.util.LinkedHashSet;

/**
 *
 * This version of AtomNeighbors matches the convention in Probe
 * where atoms bonded 1,5 are only treated as non-bonded if neither is an H
 *
 * @author mhall44
 */
public class ProbeAtomNeighbors extends AtomNeighbors {

    LinkedHashSet<Atom> neighbors15h;

    public ProbeAtomNeighbors(Atom mainAtom) {
        super(mainAtom);

        neighbors15h = new LinkedHashSet<>();

        for(Atom atom14 : neighbors14){
            for(Atom atom15 : atom14.bonds){
                if(atom15.isHydrogen() || mainAtom.isHydrogen()){
                    if( (atom15!=mainAtom) && (!neighbors12.contains(atom15))
                            && (!neighbors13.contains(atom15))
                            && (!neighbors14.contains(atom15)) ){
                        //no connection shorter than 1,5
                        neighbors15h.add(atom15);
                    }
                }
            }
        }
    }


    @Override
    public Type classifyAtom(Atom atom){
        if(neighbors15h.contains(atom))
            return Type.BONDED15H;
        else
            return super.classifyAtom(atom);
    }

}
