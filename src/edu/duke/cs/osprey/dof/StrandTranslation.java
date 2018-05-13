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

