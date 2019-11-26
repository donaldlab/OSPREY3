package edu.duke.cs.osprey.confspace.compiled;


import java.util.List;

/**
 * Manipulates atom positions in an AssignedCoords by varying continuous parameters.
 */
public interface ContinuousMotion {

	/**
	 * A description of a DegreeOfFreedom.
	 * Used to create a DegreeOfFreedom from an AssignedCoords.
	 */
	interface Description {
		ContinuousMotion build(AssignedCoords conf, ConfSpace.Pos pos);
	}

	/**
	 * Create degrees of freedom from the motion and add them to the list.
	 */
	void appendDofs(List<DegreeOfFreedom> dofs);
}
