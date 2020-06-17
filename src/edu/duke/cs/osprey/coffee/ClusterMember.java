package edu.duke.cs.osprey.coffee;

import com.hazelcast.cluster.Address;
import com.hazelcast.cluster.Member;
import com.hazelcast.config.Config;
import com.hazelcast.core.Hazelcast;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.cp.ICountDownLatch;
import com.hazelcast.instance.impl.HazelcastInstanceProxy;
import com.hazelcast.internal.serialization.impl.AbstractSerializationService;
import com.hazelcast.nio.serialization.StreamSerializer;
import com.hazelcast.spi.impl.NodeEngineImpl;
import com.hazelcast.spi.impl.operationexecutor.OperationRunner;
import com.hazelcast.spi.impl.operationservice.Operation;
import com.hazelcast.spi.impl.operationservice.impl.OperationServiceImpl;
import com.hazelcast.spi.impl.servicemanager.impl.ServiceManagerImpl;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.tools.IntRange;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.MathTools;

import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class ClusterMember implements AutoCloseable {

	/**
	 * Launch a pseudo-cluster on different threads.
	 */
	public static List<Throwable> launchPseudoCluster(int numMembers, Consumer<Cluster> block) {

		List<Throwable> exceptions = new ArrayList<>();

		var latch = new CountDownLatch(numMembers);

		var threads = IntStream.range(0, numMembers)
			.mapToObj(memberi -> {
				Thread thread = new Thread(() -> {
					try {

						// wait for all the threads to be ready
						latch.countDown();
						latch.await();

						block.accept(new Cluster("NodeDB", "job", memberi, numMembers));

					} catch (Throwable t) {
						t.printStackTrace();
						synchronized (exceptions) {
							exceptions.add(t);
						}
					}

				});
				thread.setDaemon(false);
				thread.setName("ClusterMember-" + memberi);
				thread.start();
				return thread;
			})
			.toArray(Thread[]::new);

		for (var thread : threads) {
			try {
				thread.join();
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}
		}

		return exceptions;
	}

	public final Cluster cluster;

	public final String name;
	public final HazelcastInstance inst;

	private final NodeEngineImpl nodeEngine;
	private final OperationServiceImpl opService;
	private final ServiceManagerImpl serviceManager;
	private final AbstractSerializationService serializationService;
	private final OperationRunner[] operationRunners;

	public ClusterMember(Cluster cluster) {

		this.cluster = cluster;

		name = String.format("%s-%d", cluster.name, cluster.nodeId);

		// configure the cluster
		Config cfg = new Config();
		cfg.setClusterName(cluster.id);
		cfg.setInstanceName(name);

		// disable Hazelcast's automatic phone home "feature", which is on by default
		cfg.setProperty("hazelcast.phone.home.enabled", "false");

		// let the cluster know which member is the director
		cfg.getMemberAttributeConfig().setAttribute("director", Boolean.toString(isDirector()));

		// actually start hazelcast
		inst = Hazelcast.newHazelcastInstance(cfg);

		// Hazelcast's "SPI" (service programming interface) is a total mess.
		// I think they've actually started removing it, but it's still there. Sort of.
		nodeEngine = ((HazelcastInstanceProxy)inst).getOriginal().node.getNodeEngine();
		opService = nodeEngine.getOperationService();
		serviceManager = (ServiceManagerImpl)nodeEngine.getServiceManager();
		serializationService = (AbstractSerializationService)nodeEngine.getSerializationService();
		operationRunners = nodeEngine.getOperationService().getOperationExecutor().getGenericOperationRunners();

		log("node started on cluster %s, addr=%s", cluster.id, nodeEngine.getThisAddress());
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

	private long barrierId = 0;

	public int id() {
		return cluster.nodeId;
	}

	public boolean isDirector() {
		return cluster.nodeId == 0;
	}

	public Address directorAddress() {
		return inst.getCluster().getMembers().stream()
			.filter(member -> Boolean.parseBoolean(member.getAttribute("director")))
			.findFirst()
			.orElseThrow(() -> new NoSuchElementException("can't find director cluster member"))
			.getAddress();
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

	public void barrier(long timeout, TimeUnit timeUnit) {

		ICountDownLatch latch = inst.getCPSubsystem().getCountDownLatch("barrier-" + barrierId);
		try {
			latch.trySetCount(cluster.numNodes);
			latch.countDown();
			boolean isSynced = latch.await(timeout, timeUnit);
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

	public void registerService(String name, Object service) {
		serviceManager.registerService(name, service);
	}

	public <T> void registerSerializer(Class<T> c, StreamSerializer<T> serializer) {
		serializationService.register(c, serializer);
	}

	public void sendToOthers(Supplier<Operation> op) {
		// TODO: isn't there some kind of IP broadcast we can use here?
		for (var addr : otherMemberAddresses()) {
			sendTo(op.get(), addr);
		}
	}

	public void sendTo(Operation op, Address address) {
		nodeEngine.getOperationService().invokeOnTarget(null, op, address);
		// NOTE: the future returned from invokeXXX() isn't useful unless the operation sends a response back
		// we don't really use responses, so just ignore the future
	}

	public <T> T requestFrom(Operation op, Address address, long timeout, TimeUnit timeUnit) {

		// just in case...
		if (!op.returnsResponse()) {
			throw new IllegalArgumentException(String.format("operation %s does not return a response",
				op.getClass().getSimpleName()
			));
		}

		try {
			var future = nodeEngine.getOperationService().invokeOnTarget(null, op, address);
			@SuppressWarnings("unchecked")
			T response = (T)future.get(timeout, timeUnit);
			return response;
		} catch (InterruptedException | ExecutionException | java.util.concurrent.TimeoutException ex) {
			throw new RuntimeException(String.format("%s request to %s failed",
				op.getClass().getSimpleName(), address
			), ex);
		}
	}

	public Address address() {
		return nodeEngine.getThisAddress();
	}

	public List<Address> otherMemberAddresses() {
		// TODO: cache this somehow? and update when cluster membership changes?
		return inst.getCluster().getMembers().stream()
			.filter(member -> !member.getAddress().equals(address()))
			.sorted(Comparator.comparing(Member::getUuid))
			.map(Member::getAddress)
			.collect(Collectors.toList());
	}

	public long finishedOperations() {
		return Arrays.stream(operationRunners)
			.mapToLong(runner -> runner.executedOperationsCount())
			.sum();
	}

	public void waitForOperation(long operationNum, long timeout, TimeUnit timeUnit) {

		long startNs = System.nanoTime();
		long stopNs = startNs + timeUnit.toNanos(timeout);

		while (finishedOperations() < operationNum) {

			// operation not done yet, stop waiting?
			if (System.nanoTime() >= stopNs) {
				throwTimeout("timed out waiting for operations to finish");
			}

			// keep waiting
			ThreadTools.sleep(50);
		}
	}
}
