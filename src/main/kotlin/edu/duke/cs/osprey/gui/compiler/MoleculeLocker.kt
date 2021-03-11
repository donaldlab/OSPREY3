package edu.duke.cs.osprey.gui.compiler

import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.tools.associateIdentity
import java.util.concurrent.atomic.AtomicInteger


/**
 * When the GUI is listening in on the compiler progress,
 * we end up trying to render the molecule while it's being modified by the compiler.
 *
 * We can solve that problem by making the compiler lock molecules before modifying them.
 */
class MoleculeLocker(lockers: List<Locker> = emptyList()) {

	class Locker(
		val mol: Molecule,
		/** The object on which to synchronize before using the molecule */
		val sync: Any
	) {
		/** incremented every time the molecule is unlocked */
		val sequence = AtomicInteger(0)
	}

	private val lockers =
		lockers.associateIdentity { it.mol to it }

	fun <R> lock(mol: Molecule, block: () -> R): R {

		val locker = lockers[mol]
		return if (locker != null) {

			// need to synchronize this mol
			val result = synchronized(locker.sync) {
				block()
			}
			locker.sequence.incrementAndGet()

			result

		} else {

			// no syncronization needed
			block()
		}
	}
}