package edu.duke.cs.osprey.coffee;

import com.hazelcast.cluster.Address;
import com.hazelcast.cluster.Member;
import com.hazelcast.config.Config;
import com.hazelcast.core.Hazelcast;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.cp.ICountDownLatch;
import com.hazelcast.instance.impl.HazelcastInstanceProxy;
import com.hazelcast.spi.impl.operationservice.Operation;
import com.hazelcast.spi.impl.operationservice.impl.InvocationFuture;
import com.hazelcast.spi.impl.operationservice.impl.OperationServiceImpl;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.tools.IntRange;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.MathTools;

import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;


public class ClusterMember implements AutoCloseable {

	public final Cluster cluster;

	public final String name;
	public final HazelcastInstance inst;

	private final OperationServiceImpl service;

	public ClusterMember(Cluster cluster) {

		this.cluster = cluster;

		name = String.format("%s-%d", cluster.name, cluster.nodeId);

		// configure the cluster
		Config cfg = new Config();
		cfg.setClusterName(cluster.id);
		cfg.setInstanceName(name);

		// disable Hazelcast's automatic phone home "feature", which is on by default
		cfg.setProperty("hazelcast.phone.home.enabled", "false");

		inst = Hazelcast.newHazelcastInstance(cfg);
		service = ((HazelcastInstanceProxy)inst).getOriginal().node.getNodeEngine().getOperationService();

		log("node started on cluster %s", cluster.id);
	}

	@Override
	public void close() {
		inst.getLifecycleService().shutdown();
		log("node finished");
	}

	public void log(String fmt, Object ... args) {
		Log.log(name + ": " + fmt, args);
		System.out.flush();
	}

	public void log0(String fmt, Object ... args) {
		if (cluster.nodeId == 0) {
			log(fmt, args);
		}
	}

	public static void sleep(int ms) {
		try {
			Thread.sleep(ms);
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}

	private long barrierId = 0;

	public int id() {
		return cluster.nodeId;
	}

	/**
	 * Statically partition a simple range among the cluster members into contiguous blocks, one per memeber
	 */
	public IntRange simplePartition(int size) {
		int blockSize = MathTools.divUp(size, cluster.numNodes);
		return new IntRange(
			cluster.nodeId*blockSize,
			Math.min(size, (cluster.nodeId + 1)*blockSize)
		);
	}

	public void barrier() {

		ICountDownLatch latch = inst.getCPSubsystem().getCountDownLatch("barrier-" + barrierId);
		try {
			latch.trySetCount(cluster.numNodes);
			latch.countDown();
			boolean isSynced = latch.await(1, TimeUnit.MINUTES);
			if (!isSynced) {
				throwTimeout("timed out waiting for barrier");
			}
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}

		// NOTE: don't need to destroy the latch
		// apparently it gets automatically destroyed after await() finishes

		barrierId += 1;
	}

	public void throwTimeout(String msg) {
		throw new TimeoutException(msg);
	}

	public class TimeoutException extends RuntimeException {
		public TimeoutException(String msg) {
			super(name + ": " + msg);
		}
	}

	private int nextObjectId;

	public <T> InvocationFuture<T> sendToAll(Operation op) {
		return service.invokeOnPartition(op);
	}

	public <T> InvocationFuture<T> sendTo(Operation op, Address address) {
		return service.invokeOnTarget(null, op, address);
	}

	public List<Address> memberAddresses() {
		return inst.getCluster().getMembers().stream()
			.sorted(Comparator.comparing(Member::getUuid))
			.map(Member::getAddress)
			.collect(Collectors.toList());
	}
}
