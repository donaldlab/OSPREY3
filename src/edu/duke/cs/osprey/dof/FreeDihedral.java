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

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;
import edu.duke.cs.osprey.tools.RigidBodyMotion;

/**
 *
 * @author mhall44
 */
public class FreeDihedral extends DegreeOfFreedom {
    //This is a dihedral that can move freely without any compensating motions or motions of other residues
    // e.g. a sidechain dihedral in any amino acid besides proline
    //we apply the specified angle
    
    private static final long serialVersionUID = 8005789766951509654L;
    
    Residue res;//residue we're moving
    int dihedralNum;//is this chi1, chi2 or what?
    double curVal;
    
    // temp space
    double[][] dihedralCoords;

    public FreeDihedral(Residue res, int dihedralNum) {
    	
        this.dihedralCoords = new double[4][3];
        
        this.res = res;
        this.dihedralNum = dihedralNum;
        
        // NOTE: this class is often instantiated for dihedral angles that don't exist yet
        // so don't require valid dihedrals in the constructor
        if (isValid()) {
        	this.curVal = measureDihedralDegrees();
        } else {
        	this.curVal = Double.NaN;
        }
    }
    
    public boolean isValid() {
    	return res.template.numDihedrals > 0
    		&& dihedralNum >= 0
    		&& dihedralNum < res.template.numDihedrals;
    }
    
    public void checkValid() {
    	
    	if (res.template.numDihedrals <= 0) {
    		throw new IllegalArgumentException(String.format("residue " + res.fullName + " doesn't have any dihedral angles"));
    	}
    	
    	if (dihedralNum < 0 || dihedralNum >= res.template.numDihedrals) {
    		throw new IllegalArgumentException(String.format("invalid dihedral number %d, expected in [0,%d]", dihedralNum, res.template.numDihedrals - 1));
    	}
    }
    
    public double getCurVal() {
        return curVal;
    }
    
    public double[][] updateDihedralCoords() {
    	
    	// once we examine the dihedral angle though, we require it to be valid
    	checkValid();
    	
        int[] dihAtomIndices = res.template.getDihedralDefiningAtoms(dihedralNum);
        for(int a=0; a<4; a++) {
            System.arraycopy(res.coords, 3*dihAtomIndices[a], dihedralCoords[a], 0, 3);
        }
        return dihedralCoords;
    }
    
    public double measureDihedralDegrees() {
        updateDihedralCoords();
        return Protractor.measureDihedral(dihedralCoords);
    }
    
    @Override
    public void apply(double angleDegrees) {
        
        // compute the target dihedral
        double angleRadians = Math.toRadians(angleDegrees);
        double sin = Math.sin(angleRadians);
        double cos = Math.cos(angleRadians);
        
        // measure the current dihedral
        // NOTE: measuring a dihedral requires evaluating an inverse cosine, which is slow
        // let's work with sines and cosines of dihedrals directly
        updateDihedralCoords();
        double measuredSinCos[] = Protractor.measureDihedralSinCos(dihedralCoords);
        
        // calc the dihedral rotation as a rigid body transformation relative to the current pose
        double dsin = sin*measuredSinCos[1] - cos*measuredSinCos[0];
        double dcos = cos*measuredSinCos[1] + sin*measuredSinCos[0];
        RigidBodyMotion dihRotation = new DihedralRotation(dihedralCoords[1], dihedralCoords[2], dsin, dcos);
        
        // rotate all the atoms that are moved by the dihedrals (i.e., everything beyond the third atom)
        for(int index : res.template.getDihedralRotatedAtoms(dihedralNum)) {
            dihRotation.transform(res.coords, index);
        }
        
        // store the orignal (unnormalized) value
        curVal = angleDegrees;
    }
    
    @Override
    public Residue getResidue() {
        return res;
    }
    
    public int getDihedralNumber() {
        return dihedralNum;
    }
    
    @Override
    public DegreeOfFreedom copy() {
        return new FreeDihedral(res, dihedralNum);
    }
    
    @Override
    public void setMolecule(Molecule val) {
        
        // match our residue to the one in the other molecule
        res = val.getResByPDBResNumber(res.getPDBResNumber());
    }
    
    @Override
    public DOFBlock getBlock(){
        return null;
    }

    @Override
    public String getName() {
        return "DIH"+res.getPDBResNumber()+"."+dihedralNum;
    }
}

