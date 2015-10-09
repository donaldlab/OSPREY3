/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.structure;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;

/**
 *
 * @author mhall44
 */
public class AtomNeighbors {
    //This is a classification of the neighbors of an atom as
    //directly bonded (1,2) or indirectly bonded (1,3 and 1,4)
    
    Atom mainAtom;//the atom whose neighbors is listed
    LinkedHashSet<Atom> neighbors12, neighbors13, neighbors14;
    
    
    public enum NEIGHBORTYPE {
            SELF, BONDED12, BONDED13, BONDED14, NONBONDED
    }//nonbonded here means not bonded 1,4 or closer
    
    
    
    public AtomNeighbors(Atom mainAtom){
        this.mainAtom = mainAtom;
        
        neighbors12 = new LinkedHashSet<>();
        //1,2-neighbors are just the atoms bonded to mainAtom
        neighbors12.addAll(mainAtom.bonds);
        
        
        neighbors13 = new LinkedHashSet<>();
        
        for(Atom bondedAtom : neighbors12){
            for(Atom atom13 : bondedAtom.bonds){
                if( (atom13!=mainAtom) && (!neighbors12.contains(atom13)) ){
                    //no connection shorter than 1,3
                    neighbors13.add(atom13);
                }
            }
        }
        
        neighbors14 = new LinkedHashSet<>();
        
        for(Atom atom13 : neighbors13){
            for(Atom atom14 : atom13.bonds){
                if( (atom14!=mainAtom) && (!neighbors12.contains(atom14))
                        && (!neighbors13.contains(atom14)) ){
                    //no connection shorter than 1,4
                    neighbors14.add(atom14);
                }
            }
        }
    }
    
    
    public NEIGHBORTYPE classifyAtom(Atom atom){
        //what kind of neighbor to mainAtom (if any) is atom?
        if(neighbors14.contains(atom))
            return NEIGHBORTYPE.BONDED14;
        else if(neighbors13.contains(atom))
            return NEIGHBORTYPE.BONDED13;
        else if(neighbors12.contains(atom))
            return NEIGHBORTYPE.BONDED12;
        else if(atom==mainAtom)
            return NEIGHBORTYPE.SELF;
        else
            return NEIGHBORTYPE.NONBONDED;
    }
    
    
    

    
    
    //for forcefield calculations we need to classify pairs of atoms 
    //(one in atoms1, one in atoms2) according to if they're nonbonded, 1-4 bonded, or neither.
    //internalE implies that atoms1 and atoms2 are the same, and we're calculating the internal energy
    //of their residue (or part of a residue).  In this case we need to avoid double-counting.
    //Otherwise they're assumed to be different and interacting with each other
    public static ArrayList<Atom[]> getPairs14( ArrayList<Atom> atoms1, ArrayList<Atom> atoms2, boolean internalE ){
        return getPairsByType(atoms1,atoms2,internalE,NEIGHBORTYPE.BONDED14);
    }
        
    public static ArrayList<Atom[]> getPairsNonBonded( ArrayList<Atom> atoms1, ArrayList<Atom> atoms2, boolean internalE ){
        return getPairsByType(atoms1,atoms2,internalE,NEIGHBORTYPE.NONBONDED);
    }
    
    
    public static ArrayList<Atom[]> getPairsByType( ArrayList<Atom> atoms1, 
            ArrayList<Atom> atoms2, boolean internalE, NEIGHBORTYPE type ){
        
        ArrayList<Atom[]> ans = new ArrayList<>();
        
        if(internalE){//atoms1 and atoms2 expected to be the same
            for(int atNum=0; atNum<atoms1.size(); atNum++){
                
                Atom atom1 = atoms1.get(atNum);
                AtomNeighbors neighbors = new AtomNeighbors(atom1);
                
                //To avoid double-counting, we count only pairs consisting of atom atNum and lower-numbered atoms
                for(int atNum2=0; atNum2<atNum; atNum2++){
                    Atom atom2 = atoms1.get(atNum2);
                    if( neighbors.classifyAtom(atom2) == type )
                        ans.add( new Atom[] {atom1,atom2} );
                }
            }
        }
        else {//consider all pairs
            for(Atom atom1 : atoms1){
                
                AtomNeighbors neighbors = new AtomNeighbors(atom1);
                
                for(Atom atom2 : atoms2){
                    if( neighbors.classifyAtom(atom2) == type )
                        ans.add( new Atom[] {atom1,atom2} );
                }
            }        
        }
        
        return ans;
    }
    
    
    //sometimes we have a list of selected atom pairs, and we want to separate out the ones that have a certain neighbory type
    public static ArrayList<Atom[]> getPairsByType( ArrayList<Atom[]> atomPairs, NEIGHBORTYPE type ){
        
        ArrayList<Atom[]> ans = new ArrayList<>();
        
        for(Atom[] pair : atomPairs){
            
            AtomNeighbors neighbors = new AtomNeighbors(pair[0]);
            if(neighbors.classifyAtom(pair[1]) == type){
                ans.add(pair);
            }
        }
        
        return ans;
    }
    
    //special versions
     public static ArrayList<Atom[]> getPairs14( ArrayList<Atom[]> atomPairs ){
        return getPairsByType(atomPairs, NEIGHBORTYPE.BONDED14);
    }
        
    public static ArrayList<Atom[]> getPairsNonBonded( ArrayList<Atom[]> atomPairs ){
        return getPairsByType(atomPairs, NEIGHBORTYPE.NONBONDED);
    }
        
      
    
}
