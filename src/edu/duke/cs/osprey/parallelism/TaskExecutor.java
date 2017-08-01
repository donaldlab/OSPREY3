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
