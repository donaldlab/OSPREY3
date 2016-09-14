/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import java.util.ArrayList;

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
        
        //measuring a dihedral requires evaluating an inverse cosine, which is slow
        //let's work with sines and cosines of dihedrals directly
        double curDihedralSC[] = Protractor.measureDihedralSinCos(curCoords);
        double sinDih = Math.sin(Math.PI*paramVal/180);
        double cosDih = Math.cos(Math.PI*paramVal/180);
        double sinDihChange = sinDih*curDihedralSC[1] - cosDih*curDihedralSC[0];
        double cosDihChange = cosDih*curDihedralSC[1] + sinDih*curDihedralSC[0];
        
        RigidBodyMotion dihRotation = new DihedralRotation( curCoords[1], curCoords[2], 
                sinDihChange, cosDihChange );
        
        /*double curDihedralVal = Protractor.measureDihedral(curCoords);//could go faster by not copying...hmm
        double dihedralChange = paramVal - curDihedralVal;
        //should we update curVal?
        
        RigidBodyMotion dihRotation = new DihedralRotation(curCoords[1], curCoords[2], dihedralChange);*/
        //rotate about third atom, axis = third-second atom (i.e. bond vector),
        //rotate by dihedralChange
        
        for(int index : rotatedAtomIndices)
            dihRotation.transform(res.coords, index);

        curVal = paramVal; // store the value
    }
    
    @Override
    public Residue getResidue() { return res; }
    
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
}
