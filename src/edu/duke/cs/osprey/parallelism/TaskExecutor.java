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
import edu.duke.cs.osprey.tools.HashCalculator;

import java.util.HashMap;
import java.util.Map;

public class TaskExecutor implements AutoCleanable {

	public interface Task<T> {

		T run();

		interface WithContext<T,C> extends Task<T> {

			int instanceId();
			T run(C ctx);

			default ContextId contextId() {
				return new ContextId(instanceId(), getClass());
			}

			default T run() {
				return run(null);
			}
		}
	}

	public static class ContextId {

		public final int instanceId;
		public final Class<?> taskClass;

		public ContextId(int instanceId, Class<?> taskClass) {
			this.instanceId = instanceId;
			this.taskClass = taskClass;
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(instanceId, taskClass.hashCode());
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof ContextId && equals((ContextId)other);
		}

		public boolean equals(ContextId other) {
			return this.instanceId == other.instanceId
				&& this.taskClass.equals(other.taskClass);
		}
	}

	public interface TaskListener<T> {
		void onFinished(T result);
	}

	public int getParallelism() {
		return 1;
	}

	private Map<ContextId,Object> contexts = new HashMap<>();

	public void putContext(ContextId ctxid, Object ctx) {
		contexts.put(ctxid, ctx);
	}

	public void putContext(int instanceId, Class<?> taskClass, Object ctx) {
		putContext(new ContextId(instanceId, taskClass), ctx);
	}

	public Object getContext(ContextId ctxid) {
		return contexts.get(ctxid);
	}

	public Object getContext(int instanceId, Class<?> taskClass) {
		return getContext(new ContextId(instanceId, taskClass));
	}

	public boolean isBusy() {
		return false;
	}

	public boolean isWorking() {
		return false;
	}

	public <T> void submit(Task<T> task, TaskListener<T> listener) {
		listener.onFinished(runTask(task));
	}

	protected <T> T runTask(Task<T> task) {
		if (task instanceof Task.WithContext) {
			Task.WithContext<T,Object> taskWithContext = (Task.WithContext<T,Object>)task;
			return taskWithContext.run(getContext(taskWithContext.contextId()));
		} else {
			return task.run();
		}
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
