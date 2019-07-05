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

package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.structure.Molecule;
import java.io.Serializable;

/**
 * A {@see Molecule} instance with associated continuous degrees of freedom.
 */
public class ParametricMolecule implements Serializable {
	
	/** The molecule manipulated by the degrees of freedom. */
	public final Molecule mol;
	
	/** Degrees of freedom. */
	public final List<DegreeOfFreedom> dofs;

	/** Bounds on the degrees of freedom */
	public final ObjectiveFunction.DofBounds dofBounds;
	
	/**
	 * A parametric molecule with bounds on its degrees of freedom
	 * SimpleConfSpace generates a ParametricMolecule for a particular conf (RC list),
	 * with the residue types and starting conformation implied by that conf
	 * It is helpful to bundle it with the conf's DOF bounds, especially for DEEPer and CATS
	 */
	public ParametricMolecule(Molecule mol, List<DegreeOfFreedom> dofs, ObjectiveFunction.DofBounds dofBounds) {
		this.mol = mol;
		this.dofs = dofs;
		this.dofBounds = dofBounds;
	}

	/**
	 * A parametric molecule with no degrees of freedom
	 */
	public ParametricMolecule(Molecule mol) {
		this(mol, new ArrayList<>(), new ObjectiveFunction.DofBounds(0));
	}
}
