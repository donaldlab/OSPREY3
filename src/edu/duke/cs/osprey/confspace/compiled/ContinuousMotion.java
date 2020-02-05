package edu.duke.cs.osprey.confspace.compiled;


import java.util.List;

/**
 * Manipulates atom positions in an AssignedCoords by varying continuous parameters.
 */
public interface ContinuousMotion {

	/**
	 * A description of a continuous motion affecting a conformation.
	 * Used to create ContinuousMotion instances from an AssignedCoords.
	 */
	interface ConfDescription {
		ContinuousMotion build(AssignedCoords conf, ConfSpace.Pos pos);
	}

	/**
	 * A description of a continuous motion affecting a molecule.
	 * Used to create ContinuousMotion instances from an AssignedCoords.
	 */
	interface MolDescription {
		ContinuousMotion build(AssignedCoords conf, int molInfoIndex);
	}

	/**
	 * Create degrees of freedom from the motion and add them to the list.
	 */
	void appendDofs(List<DegreeOfFreedom> dofs);
}
