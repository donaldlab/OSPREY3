/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class ResidueTemplate implements Serializable {
    //This class defines a residue type, like ARG or H2O or whatever
    //structures are composed entirely of residues. 
    //We need to identify the type of each residue whenever we want to move it it a residue-specific way
    //(e.g. to define sidechain dihedrals and rotamers),
    //or whenever we want to evaluate its energy (except by an ab initio method that takes only element 
    //types and coordinates)
    //We will also need these templates to mutate from one type to another
    
    //any of the information below can be null (i.e., not provided), but if so 
    //then attempts to call it will produce an error
    
    
    public String name;//e.g. ARG
    
    //"standard" residue
    //info from here will mostly be copied into any residue in the structure to which this template is assigned
    //(charges, standardized atom names, bonds, etc.)
    public Residue templateRes;
    
    //dihedral information
    public int numDihedrals = -1;
    public int dihedral4Atoms[][];//for each dihedral, list of 4 atoms defining it
    //these are indices in all our atom-wise arrays
    public ArrayList<ArrayList<Integer>> dihedralMovingAtoms;//list of atoms that move for each dihedral
    
    
    //rotamer information
    public int numRotamers = -1;
    double rotamericDihedrals[][];//numRotamers x numDihedrals array. 
    //For each rotamer, gives the ideal values of all dihedrals
    
    //note: this is just the standard set of ideal rotamers, we can modify this as desired
    
    public ResidueTemplate(Residue res, String name){
        //initializes only with info from a template file.  res won't even have coordinates yet
        templateRes = res;
        this.name = name;
    }
    
    //information on dihedrals
    public int[] getDihedralDefiningAtoms(int dihedralNum){
        //numbers (indices among the atoms in this template) of the 4 atoms defining the dihedral
        return dihedral4Atoms[dihedralNum];
    }
    
    public ArrayList<Integer> getDihedralRotatedAtoms(int dihedralNum){
        //the atoms that actually move (including the 4th of the dihedral-defining atoms)
        return dihedralMovingAtoms.get(dihedralNum);
    }
    
   
    
    
    void computeDihedralMovingAtoms(){
        //compute what atoms move when a dihedral is changed
        //Compute from dihedral4Atoms, and write answer in dihedralMovingAtoms
        //we assume atom 4 (from the dihedral4Atoms) moves and the others are fixed; 
        //thus the moving atoms consist of the atoms bonded to atom 3 (excluding atom 2)
        //and anything beyond them (compute recursively).  If any of these moving atoms bond back to atom 2 
        //or atom 1 or to another residue,
        //then the dihedral can't move freely, and we return an error.  
        dihedralMovingAtoms = new ArrayList<>();
        
        for(int dihedNum=0; dihedNum<numDihedrals; dihedNum++){
            
            ArrayList<Integer> movingAtoms = new ArrayList<>();//indices of atoms that move (in templateRes.atoms)
            Atom atom2 = templateRes.atoms.get(dihedral4Atoms[dihedNum][1]);
            Atom atom3 = templateRes.atoms.get(dihedral4Atoms[dihedNum][2]);
            
            for(Atom atomBeyond : atom3.bonds){
                if(atomBeyond!=atom2){//atom 2 of dihedral doesn't move
                    
                    movingAtomsDFS( atomBeyond, movingAtoms, dihedral4Atoms[dihedNum]);
                    //search recursively for other moving atoms
                }
            }
            
            dihedralMovingAtoms.add(movingAtoms);
        }
    }
    
    
    private void movingAtomsDFS(Atom atom, ArrayList<Integer> movingAtoms, int[] dihedAtoms){
        //Given an atom that moves when a dihedral is changed,
        //add it and any atoms distal to it to movingAtoms
        //dihedAtoms are the indices of the 4 atoms defining the dihedral
        if(atom==null)//atom is in another residue...this is an error
            throw new RuntimeException("ERROR: Trying to define non-free dihedral in rotamer library (bonds to another residue)");
        
        int atomIndex = templateRes.getAtomIndexByName(atom.name);
        
        if(atomIndex==dihedAtoms[0] || atomIndex==dihedAtoms[1])//there is a cycle within the residue
            throw new RuntimeException("ERROR: Trying to define non-free dihedral in rotamer library (dihedral is in a ring)");
        
        if(atomIndex == dihedAtoms[2])//dihedAtoms[2] doesn't move, though it is bonded to moving atoms
            return;
        
        //check that we haven't already counted atom as moving (since cycles are allowed within the set of moving atoms)
        if(movingAtoms.contains(atomIndex))
            return;
            
        //if we get here, atom is a new moving atom
        //add it to the list of movingAtoms, and search for atoms distal to it
        movingAtoms.add(atomIndex);
        for(Atom distalAtom : atom.bonds)
            movingAtomsDFS(distalAtom,movingAtoms,dihedAtoms);
    }
    
    
}
