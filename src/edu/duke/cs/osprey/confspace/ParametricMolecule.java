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

package edu.duke.cs.osprey.confspace;

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
}
