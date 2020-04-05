/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.TimeUnit;

import edu.duke.cs.tpie.Cleaner;
import edu.duke.cs.tpie.Cleaner.GarbageDetectable;


public class ThreadPoolTaskExecutor extends ConcurrentTaskExecutor implements GarbageDetectable {
	
	/**
	 * Controls task queue size.
	 * Set this to 0 to cause main thread to block until task thread is ready.
	 * Set this to >0 to "buffer" tasks so task threads don't have to wait on the main thread to start a task.
	 * >0 can sometimes be faster than 0, but only works if you know how many tasks you have in advance.
	 * otherwise, you can end up executing more tasks than you need.
	 * The best queue size to use is determined by the amount of work it takes to create a task vs execute it.
	 * Experiment to find the best values for your problem.
	 */
	public int queueSize = 0;
	
	private Threads threads;

	public ThreadPoolTaskExecutor() {
		threads = null;
	}
	
	public void start(int numThreads) {
		threads = new Threads(numThreads, queueSize);
		Cleaner.addCleaner(this, threads);
	}
	
	public void stop() {
		if (threads != null) {
			threads.clean();
			threads = null;
		}
	}
	
	public void stopAndWait(int timeoutMs) {
		if (threads != null) {
			threads.cleanAndWait(timeoutMs);
			threads = null;
		}
	}
	
	@Override
	public void clean() {
		stop();
	}
	
	@Override
	public int getParallelism() {
		return threads.size();
	}

	@Override
	public <T> void submit(Task<T> task, TaskListener<T> listener) {

		boolean wasAdded = false;
		while (!wasAdded) {

			checkException();

			wasAdded = threads.submit(400, TimeUnit.MILLISECONDS, () -> {
				try {

					// run the task
					T result = runTask(task);

					// send the result to the listener thread
					threads.submitToListener(() -> {
						taskSuccess(task, listener, result);
					});

				} catch (Throwable t) {
					taskFailure(task, listener, t);
				}
			});
		}

		// the task was started successfully, hooray!
		startedTask();
	}
}
