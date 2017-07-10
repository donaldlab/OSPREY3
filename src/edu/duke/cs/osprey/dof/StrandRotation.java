/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.RotationMatrix;
import edu.duke.cs.osprey.tools.VectorAlgebra;

/**
 *
 * @author mhall44
 */
public class StrandRotation extends DegreeOfFreedom {
    
    MoveableStrand strand;//residue we're moving
    int angleNum;//which of the three angles this is

    
    public StrandRotation(MoveableStrand strand, int angleNum) {
        this.strand = strand;
        this.angleNum = angleNum;
    }

    
    
    @Override
    public void apply(double paramVal) {
        
        double[] rotationCenter = VectorAlgebra.add(strand.initCenter, strand.curTrans);
        strand.curAngles[angleNum] = paramVal;
        
        //figure out how much to rotate, based on old and new rotation matrices for strand
        RotationMatrix oldMatrix = strand.curRotMatrix;
        strand.updateRotMatrixFromAngles();//new rotation (relative to starting orientation)
        RotationMatrix newMatrix = strand.curRotMatrix;
        RotationMatrix changeMatrix = newMatrix.multiply( oldMatrix.transpose() );//transpose = inverse for rotations
        
        //so now we just need to rotate the whole strand by changeMatrix about rotationCenter
        RigidBodyMotion motion = new RigidBodyMotion(rotationCenter, changeMatrix, rotationCenter);

        for(Residue res : strand.res)
            motion.transform(res.coords);
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
        return "STRROT"+strand.res.get(0).getPDBResNumber()+"."+angleNum;
    }
}