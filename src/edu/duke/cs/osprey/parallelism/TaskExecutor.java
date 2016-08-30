package edu.duke.cs.osprey.parallelism;

public class TaskExecutor {
	
	public static interface TaskListener {
		void onFinished(Runnable task);
	}
	
	public void submit(Runnable task) {
		task.run();
	}
	
	public void submit(Runnable task, TaskListener listener) {
		task.run();
		listener.onFinished(task);
	}
	
	public void waitForFinish() {
		// nothing to do
	}
}
