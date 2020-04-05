package edu.duke.cs.osprey.parallelism;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.collection.IQueue;
import com.hazelcast.config.Config;
import com.hazelcast.core.*;
import com.hazelcast.map.IMap;
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

	private static final String AliveName = "alive";
	private static final String TasksScatterName = "tasks-scatter";
	private static final String TasksGatherName = "tasks-gather";

	public final String name;
	public final int nodeId;
	public final int numNodes;
	public final Parallelism parallelism;
	public final boolean clientIsMember;

	public static final boolean DefaultClientIsMember = true;

	public Cluster(String name, int nodeId, int numNodes, Parallelism parallelism) {
		this(name, nodeId, numNodes, parallelism, DefaultClientIsMember);
	}

	public Cluster(String name, int nodeId, int numNodes, Parallelism parallelism, boolean clientIsMember) {
		this.name = name;
		this.nodeId = nodeId;
		this.numNodes = numNodes;
		this.parallelism = parallelism;
		this.clientIsMember = clientIsMember;
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

	public class Member extends ConcurrentTaskExecutor {

		private final int timeoutS;
		private final HazelcastInstance inst;

		public static final int DefaultTimeoutS = 60*5; // 5 minutes

		public Member() {
			this(DefaultTimeoutS);
		}

		public Member(int timeoutS) {

			this.timeoutS = timeoutS;

			// configure the cluster
			Config cfg = new Config();
			cfg.setClusterName(name);
			cfg.setInstanceName(String.format("%s-%d", name, nodeId));
			cfg.getQueueConfig(TasksScatterName).setMaxSize(numMembers()*2);

			// disable Hazelcast's automatic phone home "feature", which is on by default
			cfg.setProperty("hazelcast.phone.home.enabled", "false");

			inst = Hazelcast.newHazelcastInstance(cfg);
		}

		private void log(String fmt, Object ... args) {
			Log.log(inst.getName() + ": " + fmt, args);
			System.out.flush();
		}

		@Override
		public int getParallelism() {
			return parallelism.getParallelism();
		}

		@Override
		public <T> void submit(TaskExecutor.Task<T> task, TaskListener<T> listener) {
			throw new UnsupportedOperationException("Member nodes accept tasks from the Client node, not the TaskExecutor interface");
		}

		@Override
		public void putContext(ContextId ctxid, Object ctx) {
			Cluster.Task.putContext(inst, ctxid, ctx);
		}

		@Override
		public Object getContext(ContextId ctxid) {
			return Cluster.Task.getContext(inst, ctxid);
		}

		@Override
		public void clean() {
			run();
		}

		private void run() {

			log("Member node ready");

			// start a thread pool
			try (Threads threads = new Threads(parallelism.getParallelism(), 0)) {

				IQueue<Cluster.Task<?,Object>> scatter = inst.getQueue(TasksScatterName);
				IQueue<TaskResult<?>> gather = inst.getQueue(TasksGatherName);

				// wait as long as the cluster is alive
				IMap<Integer,Boolean> alive = inst.getMap(AliveName);
				while (isAlive(alive, timeoutS)) {

					// look for tasks in the queue to process
					Cluster.Task<?,Object> task;
					try {
						task = scatter.poll(1000, TimeUnit.MILLISECONDS);
					} catch (InterruptedException ex) {
						// interrupted, stop looking for new tasks
						break;
					}

					if (task != null) {

						// process the task on the thread pool
						try {

							boolean wasAdded = false;
							while (!wasAdded) {

								// NOTE: don't use ThreadPoolExecutor.submit() to send tasks, because it won't let us block.
								// access the work queue directly instead, so we can block if the thread pool isn't ready yet.
								wasAdded = threads.queue.offer(() -> {

									// run the task
									TaskResult<?> result;
									try {
										result = TaskResult.success(task.id, runTask(task));
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

								}, 400, TimeUnit.MILLISECONDS);
							}

						} catch (InterruptedException ex) {
							throw new Error(ex);
						}
					}
				}

			} catch (Throwable t) {

				log("Error on member thread");
				t.printStackTrace(System.out);

			} finally {

				// cluster isn't alive anymore (or never was), take down this cluster node
				inst.getLifecycleService().shutdown();
				log("Member node finished");
			}
		}

		private boolean isAlive(IMap<Integer,Boolean> alive, int seconds) {

			// try to check for the alive flag for a bit
			// but eventually give up and say the cluster doesn't exist
			long stopNs = System.nanoTime() + seconds*1_000_000_000L;
			do {

				Boolean isAlive = alive.get(0);
				if (isAlive == null) {

					// found nothing, wait a bit before trying again
					sleep(1000);

				} else {

					// found something
					return isAlive;
				}

			} while (System.nanoTime() < stopNs);

			// timed out waiting
			log("Timed out waiting to find cluster, shutting down member node");
			return false;
		}
	}

	public class Client extends ConcurrentTaskExecutor {

		private final Member member;
		private final HazelcastInstance inst;

		private final IMap<Integer,Boolean> alive;
		private final IQueue<Cluster.Task<?,?>> scatter;
		private final IQueue<TaskResult<?>> gather;

		private final Map<Long,TaskAndListener<?,?>> tasks = new HashMap<>();
		private final AtomicBoolean isActive = new AtomicBoolean(true);
		private final Thread listener;

		public Client() {

			// make a member first, if needed, on another thread
			if (clientIsMember) {
				member = new Member();
				Thread memberThread = new Thread(() -> member.clean());
				memberThread.setName("ClusterClientMember");
				memberThread.start();
			} else {
				member = null;
			}

			// configure the cluster
			ClientConfig cfg = new ClientConfig();
			cfg.setClusterName(name);

			// disable Hazelcast's automatic phone home "feature", which is on by default
			cfg.setProperty("hazelcast.phone.home.enabled", "false");

			// TODO: configure deserialization protection with a whitelist
			// see: https://docs.hazelcast.org/docs/latest/manual/html-single/index.html#untrusted-deserialization-protection

			inst = HazelcastClient.newHazelcastClient(cfg);

			// make the cluster "alive"
			alive = inst.getMap(AliveName);
			alive.put(0, true);

			scatter = inst.getQueue(TasksScatterName);
			gather = inst.getQueue(TasksGatherName);

			listener = new Thread(() -> {
				try {
					while (isActive.get()) {

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
			Log.log(name + "-client: " + fmt, args);
		}

		@Override
		public void putContext(ContextId ctxid, Object ctx) {
			Cluster.Task.putContext(inst, ctxid, ctx);
			if (member != null) {
				Cluster.Task.putContext(member.inst, ctxid, ctx);
			}
		}

		@Override
		public Object getContext(ContextId ctxid) {
			return Cluster.Task.getContext(inst, ctxid);
		}

		@Override
		public void clean() {

			// turn off the listener thread
			isActive.set(false);
			try {
				listener.join();
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}

			// shut down the cluster
			alive.put(0, false);
			inst.getLifecycleService().shutdown();
		}

		@Override
		public int getParallelism() {
			return inst.getCluster().getMembers().size();
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

		private static String ctxidStr(TaskExecutor.ContextId ctxid) {
			return String.format("%d_%s", ctxid.instanceId, ctxid.taskClass.getName());
		}

		static void putContext(HazelcastInstance inst, TaskExecutor.ContextId ctxid, Object ctx) {
			inst.getUserContext().put(ctxidStr(ctxid), ctx);
		}

		static Object getContext(HazelcastInstance inst, TaskExecutor.ContextId ctxid) {
			return inst.getUserContext().get(ctxidStr(ctxid));
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
