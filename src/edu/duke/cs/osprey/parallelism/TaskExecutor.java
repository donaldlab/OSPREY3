package edu.duke.cs.osprey.parallelism;

public class TaskExecutor {
	
	public static interface Task<T> {
		T run();
	}
	
	public static interface TaskListener<T> {
		void onFinished(T result);
	}
	
	public int getParallelism() {
		return 1;
	}
	
	public <T> void submit(Task<T> task, TaskListener<T> listener) {
		T result = task.run();
		listener.onFinished(result);
	}
	
	public void waitForFinish() {
		// nothing to do
	}
}
