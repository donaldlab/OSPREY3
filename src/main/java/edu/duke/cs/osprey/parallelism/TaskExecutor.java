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

import edu.duke.cs.osprey.tools.AutoCleanable;

public class TaskExecutor implements AutoCleanable {
	
	public static interface Task<T> {
		T run();
	}
	
	public static interface TaskListener<T> {
		void onFinished(T result);
	}
	
	public int getParallelism() {
		return 1;
	}

	public boolean isBusy() {
		return false;
	}

	public boolean isWorking() {
		return false;
	}

	public <T> void submit(Task<T> task, TaskListener<T> listener) {
		T result = task.run();
		listener.onFinished(result);
	}
	
	public void waitForFinish() {
		// nothing to do
	}
	
	public static class TaskException extends RuntimeException {
		
		private static final long serialVersionUID = 8523925290195831558L;
		
		public final Task<?> task;
		public final TaskListener<?> listener;
		
		public TaskException(Task<?> task, TaskListener<?> listener, Throwable cause) {
			super("A task failed, no new tasks can be submitted", cause);
			this.task = task;
			this.listener = listener;
		}
	}

	@Override
	public void clean() {
		// do nothing by default
	}
}
