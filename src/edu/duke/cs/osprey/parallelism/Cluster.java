package edu.duke.cs.osprey.parallelism;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.collection.IQueue;
import com.hazelcast.config.Config;
import com.hazelcast.core.*;
import com.hazelcast.map.IMap;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;


/**
 * Provides parallelism across a compute cluster using Hazelcast.
 */
public class Cluster {

	private static final String AliveName = "alive";
	private static final String TasksName = "tasks";
	private static final String TasksScatterName = "tasks-scatter";
	private static final String TasksGatherName = "tasks-gather";

	private static void sleep(int ms) {
		try {
			Thread.sleep(ms);
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}

	public static class Member {

		public final String clusterName;
		public final int nodeId;
		public final int maxMembers;

		private final String instanceName;

		public Member(String clusterName, int nodeId, int maxMembers) {

			this.clusterName = clusterName;
			this.nodeId = nodeId;
			this.maxMembers = maxMembers;

			instanceName = String.format("%s-%d", clusterName, nodeId);
		}

		private void log(String fmt, Object ... args) {
			edu.duke.cs.osprey.tools.Log.log(instanceName + ": " + fmt, args);
			System.out.flush();
		}

		private HazelcastInstance makeInstance() {

			// configure the cluster
			Config cfg = new Config();
			cfg.setClusterName(clusterName);
			cfg.setInstanceName(instanceName);
			cfg.getExecutorConfig(TasksName).setPoolSize(Parallelism.getMaxNumCPUs());
			cfg.getQueueConfig(TasksScatterName).setMaxSize(maxMembers*2);

			// disable Hazelcast's automatic phone home "feature", which is on by default
			cfg.setProperty("hazelcast.phone.home.enabled", "false");

			return Hazelcast.newHazelcastInstance(cfg);
		}

		public static final int DefaultTimeoutS = 60*5; // 5 minutes

		public void run() {
			run(DefaultTimeoutS);
		}

		public void run(int timeoutS) {

			HazelcastInstance inst = makeInstance();
			log("Member node ready");

			try {

				IQueue<Task<?>> scatter = inst.getQueue(TasksScatterName);
				IQueue<TaskResult<?>> gather = inst.getQueue(TasksGatherName);

				// wait as long as the cluster is alive
				IMap<Integer,Boolean> alive = inst.getMap(AliveName);
				while (isAlive(alive, timeoutS)) {

					// look for tasks in the queue to process
					Task<?> task = scatter.poll(1000, TimeUnit.MILLISECONDS);
					if (task != null) {

						// process the task
						TaskResult<?> result;
						try {
							result = new TaskResult<>(
								task.id,
								task.run(),
								null
							);
						} catch (Throwable t) {
							result = new TaskResult<>(
								task.id,
								null,
								t
							);
						}

						// send the result back
						// TODO: what if the queue is full?
						boolean wasOffered = false;
						while (!wasOffered) {
							wasOffered = gather.offer(result, 1000, TimeUnit.MILLISECONDS);
						}
					}
				}

			} catch (InterruptedException ex) {
				// fall through to cleanup
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

	public static class Client implements AutoCloseable {

		public final String clusterName;
		public final int nodeId;
		public final int maxMembers;

		private final Member member;
		private final HazelcastInstance inst;
		private final IMap<Integer,Boolean> alive;
		private final IExecutorService tasks;

		public Client(String clusterName, int maxMembers, int nodeId) {
			this(clusterName, nodeId, maxMembers, true);
		}

		public Client(String clusterName, int nodeId, int maxMembers, boolean alsoMember) {

			this.clusterName = clusterName;
			this.nodeId = nodeId;
			this.maxMembers = maxMembers;

			// make a member first, if needed
			if (alsoMember) {
				member = new Member(clusterName, maxMembers, nodeId);
				Thread memberThread = new Thread(() -> member.run());
				memberThread.setName("ClusterClientMember");
				memberThread.setDaemon(true);
				memberThread.start();
			} else {
				member = null;
			}

			// configure the cluster
			ClientConfig cfg = new ClientConfig();
			cfg.setClusterName(clusterName);

			// disable Hazelcast's automatic phone home "feature", which is on by default
			cfg.setProperty("hazelcast.phone.home.enabled", "false");

			inst = HazelcastClient.newHazelcastClient(cfg);

			// make the cluster "alive"
			alive = inst.getMap(AliveName);
			alive.put(0, true);

			// init the task executor
			tasks = inst.getExecutorService(TasksName);
		}

		private void log(String fmt, Object ... args) {
			edu.duke.cs.osprey.tools.Log.log(clusterName + "-client: " + fmt, args);
		}

		@Override
		public void close() {
			alive.put(0, false);
			tasks.destroy();
			// TODO: this causes exceptions on the Hazelcast thread pool somehow?!?!
			//   some "invocation" in the "partition" thingy is getting aborted because the client is shutting down?
			inst.getLifecycleService().shutdown();
		}

		/**
		 * Uses Hazelcast to run tasks on the cluster, with load-balancing.
		 */
		public BalancedExecutor tasks() {
			return new BalancedExecutor();
		}

		public class BalancedExecutor extends ConcurrentTaskExecutor {

			private final IQueue<Cluster.Task<?>> scatter = inst.getQueue(TasksScatterName);
			private final IQueue<TaskResult<?>> gather = inst.getQueue(TasksGatherName);

			private final Map<Long,TaskAndListener<?>> tasks = new HashMap<>();
			private final AtomicBoolean isActive = new AtomicBoolean(true);
			private final Thread listener = new Thread(() -> {

				try {
					while (isActive.get()) {

						// get the next task result, if any
						TaskResult<?> taskResult = gather.poll(400, TimeUnit.MILLISECONDS);
						if (taskResult != null) {

							// find the task for this result
							TaskAndListener<?> tal = tasks.get(taskResult.taskId);
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

			public BalancedExecutor() {
				listener.setName("ClusterClientListener");
				listener.setDaemon(true);
				listener.start();
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
			}

			@Override
			public int getParallelism() {
				return inst.getCluster().getMembers().size();
			}

			/**
			 * Submits a task to run on the cluster.
			 * The parameter type must be Serializable.
			 */
			@Override
			public <T> void submit(TaskExecutor.Task<T> task, TaskListener<T> listener) {
				try {

					// create and store the cluster task
					Cluster.Task<T> clusterTask = Cluster.Task.cast(task);
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

		/**
		 * Uses the TaskExecutor framework of Hazelcast to run tasks on the cluster.
		 * Downside: There's apparently no load-balancing. Tasks get submitted to random member nodes.
		 */
		public RandomExecutor randomTasks() {
			return new RandomExecutor();
		}

		public class RandomExecutor extends ConcurrentTaskExecutor {

			@Override
			public int getParallelism() {
				return inst.getCluster().getMembers().size();
			}

			/**
			 * Submits a task to run on the cluster.
			 * The parameter type must be Serializable.
			 */
			@Override
			public <T> void submit(TaskExecutor.Task<T> task, TaskListener<T> listener) {

				checkException();

				tasks.submit(
					Cluster.Task.cast(task),
					new ExecutionCallback<T>() {

						@Override
						public void onResponse(T response) {
							taskSuccess(task, listener, response);
						}

						@Override
						public void onFailure(Throwable t) {
							taskFailure(task, listener, t);
						}
					}
				);

				startedTask();
			}
		}
	}

	private static long nextTaskId = 0L;

	public static abstract class Task<T> implements TaskExecutor.Task<T>, Callable<T>, Serializable, HazelcastInstanceAware {

		private transient HazelcastInstance inst;

		private final long id = nextTaskId++;

		@Override
		public void setHazelcastInstance(HazelcastInstance inst) {
			this.inst = inst;
		}

		@Override
		public T call() {
			return run();
		}

		/**
		 * Makes sure the task is a cluster task.
		 */
		public static <T> Task<T> cast(TaskExecutor.Task<T> task) {
			if (!(task instanceof Task)) {
				throw new IllegalArgumentException("the task must extend Cluster.Task");
			}
			return (Task<T>)task;
		}
	}

	private static class TaskAndListener<T> {

		public final Task<T> task;
		public final TaskExecutor.TaskListener<T> listener;

		public TaskAndListener(Task<T> task, TaskExecutor.TaskListener<T> listener) {
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
	}
}
