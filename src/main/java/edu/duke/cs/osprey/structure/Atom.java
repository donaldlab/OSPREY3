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

	/**
	 * An optimized copy method
	 * (profiling shows the usual constructors are a bit slow)
	 * but doesn't copy the bonds
	 */
    public Atom copyToRes(Residue res) {
		Atom copy = new Atom();

		// copy simple properties
		copy.name = this.name;
		copy.elementType = this.elementType;
		copy.BFactor = this.BFactor;
		copy.modelAtomNumber = this.modelAtomNumber;
		copy.forceFieldType = this.forceFieldType;
		copy.type = this.type;
		copy.charge = this.charge;
		copy.elementNumber = this.elementNumber;
		copy.radius = this.radius;
		copy.mass = this.mass;

		// init the bonds
		copy.bonds = new ArrayList<>();

		// put the copy atom in the res
		copy.res = res;
		copy.indexInRes = res.atoms.size();
		res.atoms.add(copy);

		return copy;
	}

	// private constructor just for the optimized copyToRes() method,
	// so we can bypass the other slower constructors without breaking existing code
	private Atom() {}
    
    public void addBond(Atom atom2){
        //add a bond between this and atom2
        bonds.add(atom2);
        atom2.bonds.add(this);
    }
    
    
    public double[] getCoords(){
        //x,y, and z coordinates, pulled from the residue
		int i = 3*indexInRes;
        double x = res.coords[i++];
        double y = res.coords[i++];
        double z = res.coords[i];
        
        return new double[] {x,y,z};
    }

    public void setCoords(double x, double y, double z) {
    	int i = 3*indexInRes;
    	res.coords[i++] = x;
		res.coords[i++] = y;
		res.coords[i] = z;
	}
    
    public boolean isHydrogen() {
        return elementType.equalsIgnoreCase("H");
    }
    
    public boolean isCarbon() {
        return elementType.equalsIgnoreCase("C");
    }

	public boolean isOxygen() {
		return elementType.equalsIgnoreCase("O");
	}

	@Override
	public String toString() {
    	if (res != null) {
			int i = 3*indexInRes;
			return String.format("%s %s (%.3f, %.3f, %.3f)",
				res.fullName, name,
				res.coords[i], res.coords[i + 1], res.coords[i + 2]
			);
		} else {
    		return name;
		}
	}
}
