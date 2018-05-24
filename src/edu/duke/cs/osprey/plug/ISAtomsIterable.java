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
import java.util.ArrayList;
import java.util.Iterator;

/**
 *
 * Collection of intra+shell atoms that interact with an atom at1, made to be iterable over
 * 
 * @author mhall44
 */
public class ISAtomsIterable implements Iterable<Atom> {

    ProbeAtomNeighbors neighbors;
    Residue res1;
    ArrayList<Residue> shellResidues;
    
    public ISAtomsIterable(Atom at1, Residue res1, ArrayList<Residue> shellResidues){
        //iterator over atoms in the shell and in res1 that are >3 bonds away from at1
        //(note: if time-consuming should skip faraway res--same goes for pairs)
        neighbors = new ProbeAtomNeighbors(at1);
        this.res1 = res1;
        this.shellResidues = shellResidues;
    }
    
    @Override
    public Iterator<Atom> iterator() {
        return new ISAtomsIterator();
    };
    
    
    private class ISAtomsIterator implements Iterator<Atom> {
        
        int resCounter=-1, atomCounter=-1;
        boolean done = false;
        
        private ISAtomsIterator(){
            updateCounters();
        }
                    
        @Override
        public boolean hasNext() {
            return !done;
        }

        @Override
        public Atom next() {
            Residue nextRes = (resCounter==-1) ? res1 : shellResidues.get(resCounter);
            Atom nextAtom = nextRes.atoms.get(atomCounter);
            updateCounters();
            return nextAtom;
        }

        private void updateCounters(){
            Residue nextRes = (resCounter==-1) ? res1 : shellResidues.get(resCounter);
            Atom candAtom;
            do {
                atomCounter++;
                if(atomCounter==nextRes.atoms.size()){//iterate residue
                    resCounter++;
                    atomCounter=0;
                    if(resCounter==shellResidues.size()){
                        done = true;
                        break;
                    }
                    nextRes = shellResidues.get(resCounter);
                }
                candAtom = nextRes.atoms.get(atomCounter);
            } while(neighbors.classifyAtom(candAtom)!=AtomNeighbors.Type.NONBONDED);
        }
    }
}
