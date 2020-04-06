package edu.duke.cs.osprey.parallelism;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.collection.IList;
import com.hazelcast.collection.IQueue;
import com.hazelcast.config.Config;
import com.hazelcast.core.*;
import edu.duke.cs.osprey.tools.Log;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;


/**
 * Provides parallelism across a compute cluster using Hazelcast.
 */
public class Cluster {

	private static final String TasksActiveIdName = "tasks-activeContextId";
	private static final String TasksScatterName = "tasks-scatter";
	private static final String TasksGatherName = "tasks-gather";

	public static final boolean DefaultClientIsMember = true;

	// TODO: configure deserialization protection with a whitelist?
	// see: https://docs.hazelcast.org/docs/latest/manual/html-single/index.html#untrusted-deserialization-protection

	public final String name;
	public final int nodeId;
	public final int numNodes;
	public final boolean clientIsMember;

	public final String id;

	public Cluster(String name, String jobId, int nodeId, int numNodes) {
		this(name, jobId, nodeId, numNodes, DefaultClientIsMember);
	}

	public Cluster(String name, String jobId, int nodeId, int numNodes, boolean clientIsMember) {

		this.name = name;
		this.nodeId = nodeId;
		this.numNodes = numNodes;
		this.clientIsMember = clientIsMember;

		// append the job id to the cluster name,
		// so the cluster id will be unique on the local network
		this.id = String.format("%s-%s", name, jobId);
	}

	/** Read environment variables to determine cluster properties */
	public static Cluster fromSLURM(boolean clientIsMember) {
		return new Cluster(
			"SLURM",
			System.getenv("SLURM_JOB_ID"),
			Integer.parseInt(System.getenv("SLURM_PROCID")),
			Integer.parseInt(System.getenv("SLURM_NPROCS")),
			clientIsMember
		);
	}

	public int numMembers() {
		if (clientIsMember) {
			return numNodes;
		} else {
			return numNodes - 1;
		}
	}

	private static void sleep(int ms) {
		try {
			Thread.sleep(ms);
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}

	public TaskExecutor makeTaskExecutor(Parallelism parallelism) {
		if (nodeId == 0) {
			return new Client(parallelism);
		} else {
			return new Member(parallelism);
		}
	}

	static class Value<T> {

		public final String name;

		private final IList<T> list;

		public Value(HazelcastInstance inst, String name) {
			this.name = name;
			list = inst.getList(name);
		}

		public void set(T value) {
			if (list.isEmpty()) {
				list.add(value);
			} else {
				list.set(0, value);
			}
		}

		public T get() {
			if (list.isEmpty()) {
				return null;
			}
			return list.get(0);
		}
	}

	private static class NoClusterException extends RuntimeException {}

	private static class ActiveId implements Serializable {

		int id;
		boolean isActive;

		public ActiveId(int id, boolean isActive) {
			this.id = id;
			this.isActive = isActive;
		}
	}

	public class Member extends ConcurrentTaskExecutor {

		public final Parallelism parallelism;

		public static final int DefaultTimeoutS = 60*5; // 5 minutes
		public int timeoutS = DefaultTimeoutS;

		private final String name;
		private final HazelcastInstance inst;
		private final Value<ActiveId> activeId;
		private final IQueue<Cluster.Task<?,Object>> scatter;
		private final IQueue<TaskResult<?>> gather;

		private int nextContextGroupId = 0;
		private Cluster.Task<?,Object> deferredTask = null;

		public Member(Parallelism parallelism) {

			this.parallelism = parallelism;
			this.name = String.format("%s-member-%d", Cluster.this.name, nodeId);

			// configure the cluster
			Config cfg = new Config();
			cfg.setClusterName(id);
			cfg.setInstanceName(name);
			cfg.getQueueConfig(TasksScatterName).setMaxSize(numMembers()*2);

			// disable Hazelcast's automatic phone home "feature", which is on by default
			cfg.setProperty("hazelcast.phone.home.enabled", "false");

			inst = Hazelcast.newHazelcastInstance(cfg);

			log("node started on cluster %s", id);

			activeId = new Value<>(inst, TasksActiveIdName);
			scatter = inst.getQueue(TasksScatterName);
			gather = inst.getQueue(TasksGatherName);
		}

		@Override
		public void clean() {
			inst.getLifecycleService().shutdown();
			log("node finished");
		}

		private void log(String fmt, Object ... args) {
			Log.log(name + ": " + fmt, args);
			System.out.flush();
		}

		@Override
		public int getParallelism() {
			return parallelism.getParallelism();
		}

		@Override
		public <T> void submit(TaskExecutor.Task<T> task, TaskListener<T> listener) {
			// TODO: implement me?
			//  for stuff like the reference energies
			throw new UnsupportedOperationException("Member nodes accept tasks from the Client node, not the TaskExecutor interface");
		}

		private class ContextGroup extends TaskExecutor.ContextGroup {

			private final int contextGroupId = nextContextGroupId++;

			private ActiveId getActiveId() {

				// try to check for the alive flag for a bit
				// but eventually give up and say the cluster doesn't exist
				long stopNs = System.nanoTime() + timeoutS*1_000_000_000L;
				do {

					ActiveId id = activeId.get();
					if (id == null) {

						// found nothing, wait a bit before trying again
						sleep(1000);

					} else {

						// found something
						return id;
					}

				} while (System.nanoTime() < stopNs);

				// timed out waiting
				throw new NoClusterException();
			}

			@Override
			public void close() {

				// start a thread pool
				try (Threads threads = new Threads(parallelism.getParallelism(), 0)) {

					// handle any deferred tasks from previous groups
					if (deferredTask != null) {
						processTask(threads, deferredTask);
						deferredTask = null;
					}

					while (true) {

						// stop processing this group if it's inactive, or another later group has activated
						ActiveId activeId = getActiveId();
						if ((activeId.id == contextGroupId && !activeId.isActive) || activeId.id > contextGroupId) {
							break;
						}

						// look for tasks in the queue to process
						Cluster.Task<?,Object> task;
						try {
							task = scatter.poll(1000, TimeUnit.MILLISECONDS);
						} catch (InterruptedException ex) {
							// interrupted, stop looking for new tasks
							break;
						}

						// did the group change while we were waiting for this task?
						activeId = getActiveId();
						if (activeId.id > contextGroupId) {
							// yup, defer this task until we get to the next context
							deferredTask = task;
							break;
						}

						if (task != null) {
							processTask(threads, task);
						}
					}

				} catch (Throwable t) {
					throw new RuntimeException(
						String.format("%s: error processing member node context group %d",
							inst.getName(), contextGroupId
						),
						t
					);
				}
			}

			private void processTask(Threads threads, Cluster.Task<?,Object> task) {

				// process the task on the thread pool
				threads.submitLoop(400, TimeUnit.MILLISECONDS, () -> {

					// run the task
					TaskResult<?> result;
					try {
						result = TaskResult.success(task.id, task.run(getContext(task)));
					} catch (Throwable t) {
						result = TaskResult.failure(task.id, t);
					}

					// send the result back
					// TODO: what if the queue is full?
					try {
						boolean wasOffered = false;
						while (!wasOffered) {
							wasOffered = gather.offer(result, 1000, TimeUnit.MILLISECONDS);
						}
					} catch (InterruptedException ex) {
						// we're in a thread pool, just abort this task
					} catch (Exception ex) {
						log("ERROR: can't send task result to client node");
						ex.printStackTrace(System.out);
					}
				});
			}
		}

		@Override
		public ContextGroup contextGroup() {
			return new ContextGroup();
		}
	}

	public class Client extends ConcurrentTaskExecutor {

		public final Parallelism parallelism;

		private final String name;
		private final Member member;
		private final HazelcastInstance inst;

		private final Value<ActiveId> activeId;
		private final IQueue<Cluster.Task<?,?>> scatter;
		private final IQueue<TaskResult<?>> gather;

		private final Map<Long,TaskAndListener<?,?>> tasks = new HashMap<>();
		private final AtomicBoolean listenerActive = new AtomicBoolean(true);
		private final Thread listener;

		private int nextContextGroupId = 0;

		public Client(Parallelism parallelism) {

			this.parallelism = parallelism;

			this.name = String.format("%s-client-%d", Cluster.this.name, nodeId);

			// make a member first, if needed
			if (clientIsMember) {
				member = new Member(parallelism);
			} else {
				member = null;
			}

			// configure the cluster
			ClientConfig cfg = new ClientConfig();
			cfg.setClusterName(id);
			cfg.setInstanceName(name);

			// disable Hazelcast's automatic phone home "feature", which is on by default
			cfg.setProperty("hazelcast.phone.home.enabled", "false");

			log("node looking for cluster named %s ...", id);

			inst = HazelcastClient.newHazelcastClient(cfg);

			log("node started on cluster %s", id);

			activeId = new Value<>(inst, TasksActiveIdName);
			scatter = inst.getQueue(TasksScatterName);
			gather = inst.getQueue(TasksGatherName);

			listener = new Thread(() -> {
				try {
					while (listenerActive.get()) {

						// get the next task result, if any
						TaskResult<?> taskResult = gather.poll(400, TimeUnit.MILLISECONDS);
						if (taskResult != null) {

							// find the task for this result
							TaskAndListener<?,?> tal = tasks.get(taskResult.taskId);
							if (tal == null) {
								log("WARNING: received result for unknown task: %d", taskResult.taskId);
								continue;
							}

							try {

								// handle the task result
								if (taskResult.t != null) {
									taskFailure(tal.task, tal.listener, taskResult.t);
								} else {
									taskSuccessCoerceTypes(tal.task, tal.listener, taskResult.result);
								}

							} finally {
								// clean up the handled task
								tasks.remove(taskResult.taskId);
							}
						}
					}
				} catch (InterruptedException ex) {
					// exit the thread
				} catch (Throwable t) {
					t.printStackTrace(System.err);
				}
			});
			listener.setName("ClusterClientListener");
			listener.setDaemon(true);
			listener.start();
		}

		private void log(String fmt, Object ... args) {
			Log.log(name + ": " + fmt, args);
		}

		private class ContextGroup extends TaskExecutor.ContextGroup {

			private final int contextGroupId = nextContextGroupId++;

			private Cluster.Member.ContextGroup memberGroup;
			private Thread memberGroupThread;

			private ContextGroup() {

				if (member != null) {
					memberGroup = member.contextGroup();

					// start processing for the member group on another thread
					memberGroupThread = new Thread(() -> memberGroup.close());
					memberGroupThread.setName("ClusterClientMember");
					memberGroupThread.start();
				}

				// tell the cluster this context group is active now
				activeId.set(new ActiveId(contextGroupId, true));
			}

			@Override
			public void putContext(ContextId ctxid, Object ctx) {
				super.putContext(ctxid, ctx);

				// also put contexts to the member group
				if (memberGroup != null) {
					memberGroup.putContext(ctxid, ctx);
				}
			}

			@Override
			public void close() {

				// tell the cluster this group is done
				activeId.set(new ActiveId(contextGroupId, false));

				// wait for the member group to finish too
				if (memberGroupThread != null) {
					try {
						memberGroupThread.join();
					} catch (InterruptedException ex) {
						throw new RuntimeException(ex);
					}
				}
			}
		}

		@Override
		public ContextGroup contextGroup() {
			return new ContextGroup();
		}

		@Override
		public void clean() {

			// turn off the listener thread
			listenerActive.set(false);
			try {
				listener.join();
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}

			// shut down the client node
			inst.getLifecycleService().shutdown();

			// turn off the member node, if needed
			if (member != null) {
				member.clean();
			}
		}

		@Override
		public int getParallelism() {
			return numMembers();
		}

		/**
		 * Submits a task to run on the cluster.
		 * The parameter type must be Serializable.
		 * The task must extend Cluster.Task or an exception will be thrown.
		 */
		@Override
		public <T> void submit(TaskExecutor.Task<T> task, TaskListener<T> listener) {
			try {

				// create and store the cluster task
				Cluster.Task<T,Object> clusterTask = Cluster.Task.cast(task);
				tasks.put(clusterTask.id, new TaskAndListener<>(clusterTask, listener));

				// send the task to the cluster
				boolean wasAdded = false;
				while (!wasAdded) {
					checkException();
					wasAdded = scatter.offer(clusterTask, 400, TimeUnit.MILLISECONDS);
				}

				startedTask();

			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
	}

	private static long nextTaskId = 0L;

	public static abstract class Task<T,C> implements TaskExecutor.Task.WithContext<T,C>, Serializable, HazelcastInstanceAware {

		protected transient HazelcastInstance inst;

		public final long id = nextTaskId++;
		public final int instanceId;

		protected Task(int instanceId) {
			this.instanceId = instanceId;
		}

		@Override
		public void setHazelcastInstance(HazelcastInstance inst) {
			this.inst = inst;
		}

		@Override
		public int instanceId() {
			return instanceId;
		}

		/**
		 * Makes sure the task is a cluster task.
		 */
		public static <T,C> Task<T,C> cast(TaskExecutor.Task<T> task) {
			if (!(task instanceof Task)) {
				throw new IllegalArgumentException("the task must extend Cluster.Task");
			}
			return (Task<T,C>)task;
		}
	}

	private static class TaskAndListener<T,C> {

		public final Cluster.Task<T,C> task;
		public final TaskExecutor.TaskListener<T> listener;

		public TaskAndListener(Cluster.Task<T,C> task, TaskExecutor.TaskListener<T> listener) {
			this.task = task;
			this.listener = listener;
		}
	}

	private static class TaskResult<T> implements Serializable {

		long taskId;
		T result;
		Throwable t;

		private TaskResult(long taskId, T result, Throwable t) {
			this.taskId = taskId;
			this.result = result;
			this.t = t;
		}

		public static <T> TaskResult<T> success(long taskId, T result) {
			return new TaskResult<>(taskId, result, null);
		}

		public static <T> TaskResult<T> failure(long taskId, Throwable t) {
			return new TaskResult<>(taskId, null, t);
		}

		@Override
		public String toString() {
			return String.format("TaskResult[id=%d, %s]", taskId, t == null ? "success" : "failure");
		}
	}
}
