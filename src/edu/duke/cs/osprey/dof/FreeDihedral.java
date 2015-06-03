/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.RotationMatrix;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class FreeDihedral extends DegreeOfFreedom {
    //This is a dihedral that can move freely without any compensating motions or motions of other residues
    // e.g. a sidechain dihedral in any amino acid besides proline
    //we apply the specified angle
    
    Residue res;//residue we're moving
    int dihedralNum;//is this chi1, chi2 or what?
    double curVal;

    public FreeDihedral(Residue res, int dihedralNum) {
        this.res = res;
        this.dihedralNum = dihedralNum;
    }

    public double getCurVal() { return curVal; }
    
    
    @Override
    public void apply(double paramVal) {
            	
        //get indices (within res) of the 4 atoms defining this dihedral
        int[] dihAtomIndices = res.template.getDihedralDefiningAtoms(dihedralNum);
        
        //get indices of all the atoms that are moved by the dihedrals (i.e., everything beyond the third atom)
        ArrayList<Integer> rotatedAtomIndices = res.template.getDihedralRotatedAtoms(dihedralNum);
        
        
        //get current coordinates of our four atoms
        double curCoords[][] = new double[4][3];
        for(int a=0; a<4; a++)
            System.arraycopy(res.coords, 3*dihAtomIndices[a], curCoords[a], 0, 3);
        
        double curDihedralVal = Protractor.measureDihedral(curCoords);//could go faster by not copying...hmm
        
        double dihedralChange = paramVal - curDihedralVal;
        //should we update curVal?
        
        double[] dihBondVector = VectorAlgebra.subtract(curCoords[2],curCoords[1]);
        RigidBodyMotion dihRotation = new RigidBodyMotion (curCoords[2], dihBondVector, dihedralChange, false);
        //rotate about third atom, axis = third-second atom (i.e. bond vector),
        //rotate by dihedralChange
        
        for(int index : rotatedAtomIndices)
            dihRotation.transform(res.coords, index);

    	curVal = paramVal; // store the value
    }
    
    public Residue getResidue() { return res; }
    
    
    
    
}
