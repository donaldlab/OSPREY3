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

