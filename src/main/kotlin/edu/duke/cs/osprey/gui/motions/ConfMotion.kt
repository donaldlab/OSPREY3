package edu.duke.cs.osprey.gui.motions

import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.prep.DesignPosition


/**
 * Manipulates a conformation with a continuous motion.
 */
interface ConfMotion {

	/**
	 * Describes a continuous motion on a conformation.
	 */
	interface Description {
		fun copyTo(pos: DesignPosition): Description
		fun make(mol: Molecule, atomResolver: ConfLib.AtomPointer.Resolver): ConfMotion
	}
}
