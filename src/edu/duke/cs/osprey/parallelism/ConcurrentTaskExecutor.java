package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;

import static edu.duke.cs.osprey.tools.Log.log;


public abstract class ConcurrentTaskExecutor extends TaskExecutor {

	private final AtomicLong numTasksStarted = new AtomicLong(0);
	private final AtomicLong numTasksFinished = new AtomicLong(0);
	private Signal taskSignal = new Signal();
	private AtomicReference<TaskException> exception = new AtomicReference<>(null);

	@Override
	abstract public int getParallelism();

	@Override
	abstract public <T> void submit(TaskExecutor.Task<T> task, TaskListener<T> listener);

	@Override
	public boolean isBusy() {
		return getNumRunningTasks() >= getParallelism();
	}

	@Override
	public boolean isWorking() {
		return getNumRunningTasks() > 0;
	}

	@Override
	public void waitForFinish() {

		long numTasks = numTasksStarted.get();

		while (numTasksFinished.get() < numTasks) {

			// wait a bit before checking again, unless a task finishes
			taskSignal.waitForSignal(100);
		}

		// check for exceptions
		TaskException t = exception.get();
		if (t != null) {
			throw t;
		}
	}

	public long getNumRunningTasks() {
		return numTasksStarted.get() - numTasksFinished.get();
	}

	protected <T> void taskSuccess(Task<T> task, TaskListener<T> listener, T result) {

		try {

			// run the listener
			listener.onFinished(result);

		} catch (Throwable t) {
			recordException(task, listener, t);
		}

		// tell anyone waiting that we finished a task
		finishedTask();
	}

	protected void taskSuccessCoerceTypes(Task<?> task, TaskListener<?> listener, Object result) {

		try {

			// need to work around the compiler's type system here
			@SuppressWarnings("unchecked")
			TaskListener<Object> listenerObj = (TaskListener<Object>)listener;

			// run the listener
			listenerObj.onFinished(result);

		} catch (Throwable t) {
			recordException(task, listener, t);
		}

		// tell anyone waiting that we finished a task
		finishedTask();
	}

	protected void taskFailure(Task<?> task, TaskListener<?> listener, Throwable t) {
		recordException(task, listener, t);

		// the task failed, but still report finish
		finishedTask();
	}

	protected void checkException() {

		if (exception.get() != null) {
			// waitForFinish will throw the exception,
			// but after we let the rest of the pending work finish
			waitForFinish();
		}
	}

	protected void recordException(Task<?> task, TaskListener<?> listener, Throwable t) {

		// record the exception, but don't overwrite any existing exceptions
		// TODO: keep a list of all exceptions?
		exception.compareAndSet(null, new TaskException(task, listener, t));
	}

	protected void startedTask() {
		numTasksStarted.incrementAndGet();
	}

	protected void finishedTask() {
		numTasksFinished.incrementAndGet();
		taskSignal.sendSignal();
	}
}
