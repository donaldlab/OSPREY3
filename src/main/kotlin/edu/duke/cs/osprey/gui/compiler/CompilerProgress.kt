package edu.duke.cs.osprey.gui.compiler

import edu.duke.cs.osprey.tools.Stopwatch
import java.lang.Thread.sleep


class CompilerProgress(
	vararg val tasks: Task,
	private val threadGetter: () -> Thread
) {

	class Task(
		val name: String,
		val size: Int
	) {

		// NOTE: volatile is the lazy person's thread synchronization
		@Volatile var progress = 0

		fun increment() {
			progress += 1
		}

		val fraction get() =
			progress.toFloat()/size.toFloat()
	}

	@Volatile var report: ConfSpaceCompiler.Report? = null

	/**
	 * Blocks the current thread until the compiler thread finishes.
	 */
	fun waitForFinish() {
		threadGetter().join()
	}

	/**
	 * Waits for the compiler to finish,
	 * but prints progress messages to stdout while waiting
	 */
	fun printUntilFinish(intervalMs: Long = 1000) {

		val sw = Stopwatch().start()

		Thread {
			while (threadGetter().isAlive) {
				sleep(intervalMs)
				val percent = 100.0*tasks.sumBy { it.progress }/tasks.sumBy { it.size }
				println("Compiling: %.1f%%".format(percent))
			}
		}.start()

		waitForFinish()

		println("Compilation finished in ${sw.getTime(2)}")
	}
}
