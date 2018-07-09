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
