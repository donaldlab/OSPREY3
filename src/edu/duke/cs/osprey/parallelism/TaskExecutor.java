package edu.duke.cs.osprey.parallelism;

public class TaskExecutor {
	
	public static interface TaskListener {
		void onFinished(Runnable task);
	}
	
	public int getParallelism() {
		return 1;
	}
	
	public void submit(Runnable task) {
		task.run();
	}
	
	public void submit(Runnable task, TaskListener listener) {
		task.run();
		listener.onFinished(task);
	}
	
	public void waitForSpace() {
		// nothing to do
	}
	
	public void waitForFinish() {
		// nothing to do
	}
	
	public static interface NeedsCleanup {
		void cleanup();
	}
}
