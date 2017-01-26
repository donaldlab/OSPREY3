package edu.duke.cs.osprey.confspace;

import java.util.List;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.structure.Molecule;

public class ParametricMolecule {
	
	public final Molecule mol;
	public final List<DegreeOfFreedom> dofs;
	
	public ParametricMolecule(Molecule mol, List<DegreeOfFreedom> dofs) {
		this.mol = mol;
		this.dofs = dofs;
	}
}
