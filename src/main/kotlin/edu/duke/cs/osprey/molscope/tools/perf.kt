package edu.duke.cs.osprey.molscope.tools


inline fun <R> time(label: String, block: () -> R): R {
	val startNs = System.nanoTime()
	val result = block()
	val elapsedNs = System.nanoTime() - startNs
	val time = when (elapsedNs) {
		in 0 until 1_000 -> "%d ns".format(elapsedNs)
		in 0 until 1_000_000 -> "%.2f us".format(elapsedNs.toDouble()/1_000)
		in 0 until 1_000_000_000 -> "%.2f ms".format(elapsedNs.toDouble()/1_000_000)
		else -> "%.2f s".format(elapsedNs.toDouble()/1_000_000_000)
	}
	println("$label: $time")
	return result
}
