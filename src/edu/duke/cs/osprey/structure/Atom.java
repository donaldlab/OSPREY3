/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.structure;

import edu.duke.cs.osprey.tools.PeriodicTable;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class Atom implements Serializable {
    //this atom is part of a residue.  It can be assigned a forcefield type but doesn't have to be.
    //coordinates are stored by the residue, not here
    
    public Residue res;//what residue is it in?
    public int indexInRes;//what is the index of this atom in the residue?
    
    //These can be gotten from the PDB file
    public String name;//like CA or whatever
    public String elementType;//like C or whatever
    public double BFactor=0;//optional...leave all at 0 if not provided
    public int modelAtomNumber=-1;//number in PDB file
    
    
    
    //OPTIONAL (get from ResTemplate if appropriate)
    public String forceFieldType;
    public int type;//integer version
    public double charge;
    public ArrayList<Atom> bonds = new ArrayList<>();//what atoms is it bonded to?
    //the following are inferred by default from the element type (using PeriodicTable), but can be changed if needed
    public int elementNumber;
    public double radius;
    public double mass;
    

    

    
    
    public Atom(String name, String elementType) {
        this.name = name;
        PeriodicTable.setElementProperties(this,elementType);
    }
    
    public Atom(String name){
        //infer element type from name
        this.name = name;
        PeriodicTable.inferElementProperties(this);
    }
    
    
    
    public Atom copy(){
        //copy the atom, leaving the residue and bonds blank
        //(since these will need to be pointed to the right objects)
        Atom ans = new Atom(name,elementType);
        ans.BFactor = BFactor;
        ans.modelAtomNumber = modelAtomNumber;
        
        ans.indexInRes = indexInRes;
        //this is useful because we'll usually be copying from the template
        
        ans.res = null;
        
        ans.charge = charge;
        ans.type = type;
        ans.forceFieldType = forceFieldType;
        
        return ans;
    }
    
    
    public void addBond(Atom atom2){
        //add a bond between this and atom2
        bonds.add(atom2);
        atom2.bonds.add(this);
    }
    
    
    public double[] getCoords(){
        //x,y, and z coordinates, pulled from the residue
        double x = res.coords[3*indexInRes];
        double y = res.coords[3*indexInRes+1];
        double z = res.coords[3*indexInRes+2];
        
        return new double[] {x,y,z};
    }
    
    public boolean isHydrogen() {
        return elementType.equalsIgnoreCase("H");
    }
    
    public boolean isCarbon() {
        return elementType.equalsIgnoreCase("C");
    }

}