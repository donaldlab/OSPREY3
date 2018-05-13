/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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
