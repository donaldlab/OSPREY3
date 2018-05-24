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

package edu.duke.cs.osprey.plug;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.ProbeAtomNeighbors;
import edu.duke.cs.osprey.structure.Residue;
import java.util.Iterator;

/**
 *
 * @author mhall44
 */
public class Res2AtomsIterable implements Iterable<Atom> {
    
    ProbeAtomNeighbors neighbors;
    Residue res1;
    Residue res2;
    
    public Res2AtomsIterable(Atom at1, Residue res1, Residue res2){
        //iterator over atoms in res2 that are >3 bonds away from at1
        neighbors = new ProbeAtomNeighbors(at1);//DEBUG!!  really only need this when at1, res2 close enough
        this.res1 = res1;
        this.res2 = res2;
    }
    
    @Override
    public Iterator<Atom> iterator() {
        return new Res2AtomsIterator();
    };
        
    private class Res2AtomsIterator implements Iterator<Atom> {
        
        int atomCounter=-1;
        boolean done = false;
        
        private Res2AtomsIterator(){
            updateCounter();
        }
        
                    
        @Override
        public boolean hasNext() {
            return !done;
        }

        @Override
        public Atom next() {
            Atom nextAtom = res2.atoms.get(atomCounter);
            updateCounter();
            return nextAtom;
        }

        private void updateCounter(){
            //update counters
            Atom candAtom;
            do {
                atomCounter++;
                if(atomCounter==res2.atoms.size()){//iterate residue
                    done = true;
                    break;
                }
                candAtom = res2.atoms.get(atomCounter);
            } while(neighbors.classifyAtom(candAtom)!=AtomNeighbors.Type.NONBONDED);
        }
    }
}
