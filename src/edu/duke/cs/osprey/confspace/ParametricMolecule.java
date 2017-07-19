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
