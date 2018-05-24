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

package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Residue;

/**
 *
 * @author mhall44
 */
public class StrandTranslation extends DegreeOfFreedom {
    
    
    MoveableStrand strand;//residue we're moving
    int coordNum;//Is this translation along the x, y, or z axis?

    
    public StrandTranslation(MoveableStrand strand, int coordNum) {
        this.strand = strand;
        this.coordNum = coordNum;
    }

    
    
    @Override
    public void apply(double paramVal) {
        
        double motion = paramVal - strand.curTrans[coordNum];//how far we need to move
        
        for(Residue res : strand.res){
            int numAtoms = res.atoms.size();
            for(int atomNum=0; atomNum<numAtoms; atomNum++){
                res.coords[3*atomNum+coordNum] += motion;
            }
        }
        
        strand.curTrans[coordNum] = paramVal;
    }
    
    public MoveableStrand getMoveableStrand(){
        return strand;
    }
    
    @Override
    public DOFBlock getBlock(){
        return strand;
    }
    
    @Override
    public String getName() {
        return "STRTRANS"+strand.res.get(0).getPDBResNumber()+"."+coordNum;
    }
}
