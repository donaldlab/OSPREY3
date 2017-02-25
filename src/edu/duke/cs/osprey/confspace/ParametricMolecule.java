package edu.duke.cs.osprey.confspace;

import java.util.List;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.structure.Molecule;

/**
 * A {@see Molecule} instance with associated continuous degrees of freedom.
 */
public class ParametricMolecule {
	
	/** The molecule manipulated by the degrees of freedom. */
	public final Molecule mol;
	
	/** Degrees of freedom. */
	public final List<DegreeOfFreedom> dofs;
	
	/** */
	public ParametricMolecule(Molecule mol, List<DegreeOfFreedom> dofs) {
		this.mol = mol;
		this.dofs = dofs;
	}
}
