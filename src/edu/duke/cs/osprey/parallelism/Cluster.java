package edu.duke.cs.osprey.parallelism;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.config.Config;
import com.hazelcast.core.*;
import com.hazelcast.map.IMap;

import java.io.Serializable;
import java.util.concurrent.Callable;


/**
 * Provides parallelism across a compute cluster using Hazelcast.
 */
public class Cluster {

	private static final String AliveName = "alive";
	private static final String TasksName = "tasks";

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

		private final String instanceName;

		public Member(String clusterName, int nodeId) {

			this.clusterName = clusterName;
			this.nodeId = nodeId;

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
			cfg.getExecutorConfig(TasksName).setPoolSize(Parallelism.getMaxNumCPUs());
			cfg.setInstanceName(instanceName);

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

			// wait as long as the cluster is alive
			IMap<Integer,Boolean> alive = inst.getMap(AliveName);
			while (isAlive(alive, timeoutS)) {
				sleep(1000);
				// TODO: process a task queue while we wait?
			}

			// cluster isn't alive anymore (or never was), take down this cluster node
			inst.getLifecycleService().shutdown();
			log("Member node finished");
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

		private final HazelcastInstance inst;
		private final HazelcastInstance memberInst;
		private final IMap<Integer,Boolean> alive;
		private final IExecutorService tasks;

		public Client(String clusterName, int nodeId) {
			this(clusterName, nodeId, true);
		}

		public Client(String clusterName, int nodeId, boolean alsoMember) {

			this.clusterName = clusterName;
			this.nodeId = nodeId;

			// make a member instance first, if needed
			if (alsoMember) {
				memberInst = new Member(clusterName, nodeId).makeInstance();
			} else {
				memberInst = null;
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
			inst.getLifecycleService().shutdown();
			if (memberInst != null) {
				memberInst.getLifecycleService().shutdown();
			}
		}

		/* TODO: make a load-balanced, queue-based task executor

		public BalancedExecutor tasks() {
			return new BalancedExecutor();
		}

		public class BalancedExecutor extends ConcurrentTaskExecutor {
			// TODO
		}
		*/

		public RandomExecutor randomTasks() {
			return new RandomExecutor();
		}

		/**
		 * Uses the TaskExecutor framework of Hazelcast to run tasks on the cluster.
		 * Downside: There's apparently no load-balancing. Tasks get submitted to random member nodes.
		 */
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

				// make sure the task is a cluster task
				if (!(task instanceof Cluster.Task)) {
					throw new IllegalArgumentException("the task must extend Cluster.Task");
				}
				@SuppressWarnings("unchecked")
				Cluster.Task<T> clusterTask = (Cluster.Task<T>)task;

				checkException();

				tasks.submit(
					clusterTask,
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

	public static abstract class Task<T> implements TaskExecutor.Task<T>, Callable<T>, Serializable, HazelcastInstanceAware {

		private transient HazelcastInstance inst;

		@Override
		public void setHazelcastInstance(HazelcastInstance inst) {
			this.inst = inst;
		}

		@Override
		public T call() {
			return run();
		}
	}
}
