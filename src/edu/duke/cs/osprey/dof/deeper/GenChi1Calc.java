/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof.deeper;

import edu.duke.cs.osprey.dof.DihedralRotation;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

/**
 *
 * This class handles the calculation of the generalized chi1
 * When we're moving around sidechains as a rigid body, idealization tells us where to place the CB
 * but this still leaves one degree of freedom (chi1, when defined; else a generalized chi1)
 * to determine the sidechain pose
 * 
 * @author mhall44
 */
public class GenChi1Calc {
        
    public static double getGenChi1(Residue res){//Get the generalized chi1 of a residue
        
        if( res.template.name.equalsIgnoreCase("GLY") )
            return 0;
        else{
            String lastAtom = getGenChi1LastAtom(res.template.name);
            
            double lastCoords[] = res.getCoordsByAtomName(lastAtom);
            double NCoords[] = res.getCoordsByAtomName("N");
            double CACoords[] = res.getCoordsByAtomName("CA");
            double CBCoords[] = res.getCoordsByAtomName("CB");
            
            if(lastCoords==null || NCoords==null || CACoords==null | CBCoords==null){
                //not a protein residue.  Doesn't have gen chi1
                return Double.NaN;
            }
            
            
            //coordinates defining dihedral
            double dihCoords[][] = new double[][] {NCoords,CACoords,CBCoords,lastCoords};
            
            return Protractor.measureDihedral(dihCoords);
        }
    }
    
    
    
    public static void setGenChi1(Residue res, double ang){
        //set the residue to have the specified gen chi1 (in degrees)

        if( (!HardCodedResidueInfo.hasAminoAcidBB(res)) 
                || res.template.name.equalsIgnoreCase("GLY")
                || res.template.name.equalsIgnoreCase("PRO") ){
            //Glycine doesn't have a generalized chi1, 
            //and we cannot freely change proline's (sidechain idealization sets
            //a chi1 for proline that should remain in place)
            return;
        }

        
        double curGenChi1 = getGenChi1(res);
        double genChi1Change = ang - curGenChi1;
        
        DihedralRotation dihRot = new DihedralRotation( res.getCoordsByAtomName("CA"),
                res.getCoordsByAtomName("CB"), genChi1Change );
        
        //now carry out the rotation on all the sidechain atoms
        //(expect CB, since it's part of the bond begin rotated around)
        int numAtoms = res.atoms.size();
        for(int atomNum=0; atomNum<numAtoms; atomNum++){
            
            String atomName = res.atoms.get(atomNum).name;
            
            if(SidechainIdealizer.isSidechainAtom(atomName, res)){
                if( ! ( atomName.equalsIgnoreCase("CB") ) )
                    dihRot.transform(res.coords, atomNum);
            }
        }

    }
    
    
    
    public static String getGenChi1LastAtom(String resName){//Get atom name X where generalized chi1 is N-CA-CB-X

            if( resName.equalsIgnoreCase("val") || resName.equalsIgnoreCase("ile") )
                return "CG1";
            else if( resName.equalsIgnoreCase("ser") )
                return "OG";
            else if( resName.equalsIgnoreCase("thr") )
                return "OG1";
            else if( resName.equalsIgnoreCase("cys") || resName.equalsIgnoreCase("cyx") )
                return "SG";
            else if( resName.equalsIgnoreCase("ala") )
                return "HB1";
            else
                return "CG";
    }
    
}
