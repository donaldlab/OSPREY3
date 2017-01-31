package edu.duke.cs.osprey.parallelism;

public class TaskExecutor {
	
	public static interface TaskListener<T extends Runnable> {
		void onFinished(T task);
	}
	
	public int getParallelism() {
		return 1;
	}
	
	public void submit(Runnable task) {
		task.run();
	}
	
	public <T extends Runnable> void submit(T task, TaskListener<T> listener) {
		task.run();
		if (listener != null) {
			listener.onFinished(task);
		}
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
