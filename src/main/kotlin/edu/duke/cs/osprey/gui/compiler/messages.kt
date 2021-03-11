package edu.duke.cs.osprey.gui.compiler

import cuchaz.kludge.tools.x
import cuchaz.kludge.tools.y
import cuchaz.kludge.tools.z
import org.joml.Vector3dc


fun Throwable.msgOrName(): String = message ?: javaClass.simpleName

/**
 * Flatten the cause chain into a list.
 */
fun Throwable.causes(): List<Throwable> =
	ArrayList<Throwable>().apply {
		var cause = cause
		while (cause != null) {
			add(cause)
			cause = cause.cause
		}
	}

data class CompilerWarning(val msg: String, val extraInfo: String? = null)

class CompilerError(
	/** A short message to the user to signal what went wrong */
	val msg: String,
	/** A longer message describing more details about the error */
	val extraInfo: String? = null,
	/** Exception(s) that caused the error */
	cause: Throwable? = null
) : RuntimeException(msg, cause)

fun <T:Vector3dc> T.checkForErrors(source: String): T = apply {
	// check that all the coords are valid (ie, not NaN)
	if (!x.isFinite() || !y.isFinite() || !z.isFinite()) {
		throw CompilerError("Bad coordinates at $source: [%12.6f,%12.6f,%12.6f]".format(x, y, z))
	}
}
