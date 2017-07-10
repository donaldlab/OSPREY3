/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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