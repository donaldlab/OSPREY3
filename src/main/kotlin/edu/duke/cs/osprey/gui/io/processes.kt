package edu.duke.cs.osprey.gui.io

import java.util.concurrent.ConcurrentLinkedQueue


/**
 * Reads stdout and stderr from a process and sends the lines
 * back in a ConcurrentLinkedQueue instance.
 */
class ProcessStreamer(processBuilder: ProcessBuilder) {

	init {
		// combine stdout and stderr so we only have to read one stream
		processBuilder.redirectErrorStream(true)
	}

	// start the process
	private val process = processBuilder.start()

	// read the result in a separate thread
	val console = ConcurrentLinkedQueue<String>()

	private val thread = Thread {
		process.inputStream.bufferedReader().forEachLine { line ->
			console.add(line)
		}
	}.apply {
		name = "ProcessStreamer"
		isDaemon = true
		start()
	}

	// TODO: allow timeout?
	fun waitFor() = apply {
		process.waitFor()
		thread.join()
	}

	val exitCode get() = process.exitValue()
}

fun ProcessBuilder.stream() = ProcessStreamer(this)
