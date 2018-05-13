/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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

