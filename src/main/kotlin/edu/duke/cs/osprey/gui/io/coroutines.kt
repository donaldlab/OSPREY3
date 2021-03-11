package edu.duke.cs.osprey.gui.io

import kotlinx.coroutines.*
import kotlinx.coroutines.sync.Semaphore
import kotlin.coroutines.coroutineContext


class LaunchLimits(val maxInFlight: Int) {

	val semaphore = Semaphore(maxInFlight)

	suspend fun launch(block: suspend CoroutineScope.() -> Unit): Job {

		// use the coroutine scope of the caller
		with(CoroutineScope(coroutineContext)) {

			// use a semaphore to make sure we don't have too many coroutine launches in-flight at once
			semaphore.acquire()

			return launch {
				try {
					block()
				} finally {
					semaphore.release()
				}
			}
		}
	}
}
