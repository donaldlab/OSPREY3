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
