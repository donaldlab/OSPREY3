/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RotationMatrix;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;

/**
 *
 * This is a subportion of the molecule that can move as a rigid body
 * there may be other (e.g., dihedral) motions on top of that
 * 
 * 
 * @author mhall44
 */
public class MoveableStrand implements Serializable, DOFBlock {
    
    ArrayList<Residue> res;//list of residues in the molecule that make up this "strand"
    
    //the strand has three translation and three rotation degrees of freedom
    StrandTranslation[] translations;
    StrandRotation[] rotations;
    
    //rotations will be performed about the "center" of the strand
    //which is the initCenter plus any translations that have been performed
    double initCenter[]; 
    
    //the current state of the strand is that is has been translated by curTrans, then 
    //rotated by curRotMatrix about (initCenter+curTrans)
    double curTrans[];//current value of the three translations
    RotationMatrix curRotMatrix;
    //curRotMatrix is in turn defined by curAngles
    double curAngles[];//in degrees
    
    
    public static final double maxStrandRot = 5;//maximum strand rotation Tait-Bryan angle, in degrees
    public static final double maxStrandTrans = 1.2;//maximum strand translation in any dimension, in angstroms
    //Each of these motions can go in either direction
    
    
    public MoveableStrand(){
        
    }
    
    
    public MoveableStrand(ArrayList<Residue> res){
        //Set up a moveable strand consisting of these res.  
        //Set it up to be untranslated and unrotated in the current conformation
        this.res = res;
        
        translations = new StrandTranslation[3];
        rotations = new StrandRotation[3];
        
        for(int coord=0; coord<3; coord++){
            translations[coord] = new StrandTranslation(this,coord);
            rotations[coord] = new StrandRotation(this,coord);
        }
        
        initCenter = strandCOM();
        curTrans = new double[3];//no translation yet
        curRotMatrix = RotationMatrix.identity();//no rotation either
        curAngles = new double[3];
    }
    
    
    void updateRotMatrixFromAngles(){
        //update the curRotMatrix to reflect the curAngles
        
        //this form of the rotation matrix is designed to allow any rotation,
        //while using a 3-D parameterization (quaternions would introduce a messy equality constraint)
        //and to rotate in roughly orthogonal directions (roll, pitch, and yaw) for small rotations
        //for larger ones, there could be gimbal lock, 
        //but this is OK because we can still reach any rotation we want
        
        
        //first, rotate by the first angle about the x-axis
        RotationMatrix rot1 = new RotationMatrix(1, 0, 0, curAngles[0], false);
        //then y and z
        RotationMatrix rot2 = new RotationMatrix(0, 1, 0, curAngles[1], false);
        RotationMatrix rot3 = new RotationMatrix(0, 0, 1, curAngles[2], false);
        
        //perform rot1, then 2 and 3
        curRotMatrix = rot3.multiply( rot2.multiply(rot1) );
    }
    
    
    private double[] strandCOM(){
        //center of mass of the strand
        double com[] = new double[3];
        int numAtoms = 0;
        
        for(Residue curRes : res){
            int resNumAtoms = curRes.atoms.size();
            numAtoms += resNumAtoms;
            
            for(int atomNum=0; atomNum<resNumAtoms; atomNum++){
                for(int a=0; a<3; a++)
                    com[a] += curRes.coords[3*atomNum+a];
            }
        }
        
        for(int a=0; a<3; a++)
            com[a] /= numAtoms;
        
        return com;
    }
    
    
    public ArrayList<DegreeOfFreedom> getDOFs(){
        //return a list of this strand's 6 DOFs
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        ans.addAll(Arrays.asList(rotations));
        ans.addAll(Arrays.asList(translations));
        
        return ans;
    }
    
    public ArrayList<Residue> getResidues(){
        return res;
    }
    
    
    public static double[] getStrandDOFBounds(DegreeOfFreedom strandDOF){
        //What are the bounds on this strand rigid-motion DOF?
        if(strandDOF instanceof StrandRotation){
            return new double[] {-maxStrandRot,maxStrandRot};
        }
        else if(strandDOF instanceof StrandTranslation){
            return new double[] {-maxStrandTrans,maxStrandTrans};
        }
        else
            throw new RuntimeException("ERROR: Not a strand DOF");
    }

    @Override
    public DOFBlock copyForNewMolecule(Molecule mol, LinkedHashMap<DegreeOfFreedom, DegreeOfFreedom> copiedDOFMap) {
        
        MoveableStrand copiedStrand = new MoveableStrand();
        
        copiedStrand.res = new ArrayList<>();
        for(Residue curRes : res)
            copiedStrand.res.add( mol.getResByPDBResNumber(curRes.getPDBResNumber()) );
        
        copiedStrand.translations = new StrandTranslation[3];
        copiedStrand.rotations = new StrandRotation[3];
        for(int coord=0; coord<3; coord++){
            copiedStrand.translations[coord] = new StrandTranslation(copiedStrand,coord);
            copiedStrand.rotations[coord] = new StrandRotation(copiedStrand,coord);
            copiedDOFMap.put(translations[coord], copiedStrand.translations[coord]);
            copiedDOFMap.put(rotations[coord], copiedStrand.rotations[coord]);
            copiedStrand.translations[coord].curVal = translations[coord].curVal;
            copiedStrand.rotations[coord].curVal = rotations[coord].curVal;
        }
            
        copiedStrand.initCenter = initCenter;
        copiedStrand.curTrans = curTrans.clone();
        copiedStrand.curRotMatrix = curRotMatrix.copy();
        copiedStrand.curAngles = curAngles.clone();
        
        return copiedStrand;
    }
    
}
