package edu.duke.cs.osprey.coffee.seqdb;

import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.coffee.Serializers;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.*;
import org.mapdb.*;

import java.io.File;
import java.math.MathContext;
import java.util.*;


public class SeqDB implements AutoCloseable {

	public static final String ServiceName = "seqdb";

	public final MultiStateConfSpace confSpace;
	public final MathContext mathContext;
	public final File file;
	public final ClusterMember member;

	final DB db;
	final HTreeMap<int[],SeqInfo> sequencedSums;
	final HTreeMap<Integer,StateZ> unsequencedSums;

	public SeqDB(MultiStateConfSpace confSpace, MathContext mathContext, File file, ClusterMember member) {

		this.confSpace = confSpace;
		this.mathContext = mathContext;
		this.file = file;
		this.member = member;

		// open the DB, but only on the driver member
		if (member.isDirector()) {

			if (file != null) {
				db = DBMaker.fileDB(file)
					.fileMmapEnableIfSupported() // use memory-mapped files if possible (can be much faster)
					.make();
			} else {
				db = DBMaker.memoryDB()
					.make();
			}

			// open the tables

			sequencedSums = db.hashMap("sequenced-sums")
				.keySerializer(Serializers.mapdbSeq(confSpace.seqSpace))
				.valueSerializer(Serializers.mapdbSeqInfo(confSpace))
				.createOrOpen();

			unsequencedSums = db.hashMap("unsequenced-sums")
				.keySerializer(Serializer.INTEGER)
				.valueSerializer(Serializers.mapdbStateZ)
				.createOrOpen();

		} else {
			db = null;
			sequencedSums = null;
			unsequencedSums = null;
		}

		// register NodeDB with hazelcast
		member.registerService(ServiceName, this);
	}

	public BigMath bigMath() {
		return new BigMath(mathContext);
	}

	public Batch batch() {
		return new Batch(this);
	}

	void commitBatch(SaveOperation op) {
		// TODO: NEXTTIME: fix the concurrency issues, not sure what's going wrong yet
		synchronized (this) {

			for (var sum : op.sequencedSums) {

				// convert the sum to a SeqInfo
				SeqInfo seqInfo = new SeqInfo(confSpace);
				System.arraycopy(sum.statezs, 0, seqInfo.statezs, 0, confSpace.sequencedStates.size());

				// combine with the old sums if needed
				SeqInfo oldSeqInfo = sequencedSums.get(sum.seq);
				if (oldSeqInfo != null) {
					for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
						var statez = sum.statezs[state.sequencedIndex];
						combineSums(statez, oldSeqInfo.statezs[state.sequencedIndex]);
						fixRoundoffError(statez);
					}
				}

				sequencedSums.put(sum.seq, seqInfo);
			}

			for (var sum : op.unsequencedSums) {

				var statez = sum.statez;

				// combine with the old sum if needed
				var statezOld = unsequencedSums.get(sum.unsequencedIndex);
				if (statezOld != null) {
					combineSums(statez, statezOld);
					fixRoundoffError(statez);
				}

				unsequencedSums.put(sum.unsequencedIndex, statez);
			}

			db.commit();
		}
	}

	private void combineSums(StateZ statez, StateZ statezOld) {
		statez.zSumBounds.upper = bigMath()
			.set(statez.zSumBounds.upper)
			.add(statezOld.zSumBounds.upper)
			.get();
		statez.zSumBounds.lower = bigMath()
			.set(statez.zSumBounds.lower)
			.add(statezOld.zSumBounds.lower)
			.get();
		statez.zSumDropped = bigMath()
			.set(statez.zSumDropped)
			.add(statezOld.zSumDropped)
			.get();
	}

	private void fixRoundoffError(StateZ statez) {
		// trust the lower bound more, since it's based on minimizations
		if (!statez.zSumBounds.isValid()) {
			// TODO: throw an Exception if the error is bigger than what we'd expect from roundoff?
			statez.zSumBounds.upper = statez.zSumBounds.lower;
		}
	}

	@Override
	public void close() {
		if (db != null) {
			db.close();
		}
	}

	private void checkDirector() {
		if (!member.isDirector()) {
			throw new IllegalStateException("Tried to access SeqDB directly on non-director cluster member " + member.name);
		}
	}

	/**
	 * Returns accumulated Z values for the unsequenced state,
	 * or null if that state has not yet been explored.
	 * (you probably don't want this unless you're debugging)
	 */
	StateZ getSum(MultiStateConfSpace.State state) {

		checkDirector();

		return unsequencedSums.get(state.unsequencedIndex);
	}

	/**
	 * Returns the current Z bounds and dropped Z for the unsequenced state.
	 */
	public StateZ getUnsequenced(MultiStateConfSpace.State state) {

		checkDirector();

		// return the bound if it exists
		StateZ z = getSum(state);
		if (z != null) {
			return z;
		}

		// otherwise, return unknown
		return StateZ.makeUnknown();
	}

	/**
	 * Returns accumulated sequenced state Z values for the queried sequence,
	 * or null if the sequence has not yet been explored.
	 * (you probably don't want this unless you're debugging)
	 */
	SeqInfo getSums(Sequence seq) {

		checkDirector();

		return sequencedSums.get(seq.rtIndices);
	}

	/**
	 * Returns accumulated sequenced state Z values for all sequences
	 * (you probably don't want this unless you're debugging)
	 */
	Iterable<Map.Entry<Sequence,SeqInfo>> getSums() {

		checkDirector();

		return () -> new Iterator<>() {

			final Iterator<Map.Entry<int[],SeqInfo>> iter = sequencedSums.getEntries().iterator();

			@Override
			public boolean hasNext() {
				return iter.hasNext();
			}

			@Override
			public Map.Entry<Sequence, SeqInfo> next() {

				Map.Entry<int[],SeqInfo> entry = iter.next();

				Sequence seq = new Sequence(confSpace.seqSpace, entry.getKey());
				SeqInfo seqInfo = entry.getValue();

				return new AbstractMap.SimpleEntry<>(seq, seqInfo);
			}
		};
	}

	/**
	 * Returns the current state Z bounds for unexplored subtrees of the partially-assigned sequence.
	 * As more subtrees get explored, those Z values will be transfered to more fully-assigned sequences.
	 * Bounds for fully-explored partial sequences will be zero.
	 */
	public SeqInfo getUnexplored(Sequence seq) {

		checkDirector();

		if (seq.isFullyAssigned()) {
			throw new IllegalArgumentException("need partially-assigned sequence");
		}

		// return the bound, if it exists
		SeqInfo seqInfo = sequencedSums.get(seq.rtIndices);
		if (seqInfo != null) {
			return seqInfo;
		}

		// otherwise, make an unknown bound
		return SeqInfo.makeUnknown(confSpace);
	}

	/**
	 * Returns the current sequenced state Z bounds for the fully-assigned sequence.
	 */
	public SeqInfo get(Sequence seq) {

		checkDirector();

		if (!seq.isFullyAssigned()) {
			throw new IllegalArgumentException("need fully-assigned sequence");
		}

		// return the bound, if it exists
		SeqInfo seqInfo = getSums(seq);
		if (seqInfo != null) {
			addZAncestry(seq, seqInfo);
			return seqInfo;
		}

		// otherwise, make an unknown bound
		return SeqInfo.makeUnknown(confSpace);
	}

	/**
	 * returns Z bounds for all sequences
	 * returns bounds for both full and partial sequences
	 */
	public Iterable<Map.Entry<Sequence,SeqInfo>> getSequenced() {

		checkDirector();

		return () -> new Iterator<>() {

			final Iterator<Map.Entry<int[],SeqInfo>> iter = sequencedSums.getEntries().iterator();

			@Override
			public boolean hasNext() {
				return iter.hasNext();
			}

			@Override
			public Map.Entry<Sequence,SeqInfo> next() {

				Map.Entry<int[],SeqInfo> entry = iter.next();

				Sequence seq = new Sequence(confSpace.seqSpace, entry.getKey());
				SeqInfo seqInfo = entry.getValue();

				if (seq.isFullyAssigned()) {
					addZAncestry(seq, seqInfo);
				}

				return new AbstractMap.SimpleEntry<>(seq, seqInfo);
			}
		};
	}

	private void addZAncestry(Sequence seq, SeqInfo seqInfo) {

		int[] rtIndices = new int[confSpace.seqSpace.positions.size()];
		System.arraycopy(seq.rtIndices, 0, rtIndices, 0, rtIndices.length);

		// add uncertainty from partial sequence ancestry
		// NOTE: assumes tree pos order follows seq pos order
		for (int i=rtIndices.length - 1; i>=0; i--) {
			rtIndices[i] = Sequence.Unassigned;

			SeqInfo parentSeqInfo = sequencedSums.get(rtIndices);
			if (parentSeqInfo != null) {
				// couldn't that unexplored subtree contain no confs for this seq?
				// NOTE: don't add the lower bounds, the subtree need not necessarily contain confs for this sequence
				for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
					seqInfo.statezs[state.sequencedIndex].zSumBounds.upper = bigMath()
						.set(seqInfo.statezs[state.sequencedIndex].zSumBounds.upper)
						.add(parentSeqInfo.statezs[state.sequencedIndex].zSumBounds.upper)
						.get();
				}
			}
		}
	}

	public String dump() {

		checkDirector();

		StringBuilder buf = new StringBuilder();
		buf.append("Unsequenced");
		for (MultiStateConfSpace.State state : confSpace.unsequencedStates) {
			var statez = getUnsequenced(state);
			buf.append(String.format("\n%10s  gBounds=%s  dropped=%s",
				state.name,
				Log.formatBigEngineering(statez.zSumBounds),
				Log.formatBigEngineering(statez.zSumDropped)
			));
		}
		buf.append("\nSequenced");
		for (Map.Entry<Sequence,SeqInfo> entry : getSequenced()) {
			Sequence seq = entry.getKey();
			SeqInfo seqInfo = entry.getValue();
			buf.append(String.format("\n[%s]", seq));
			for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
				var statez = seqInfo.get(state);
				buf.append(String.format("\n%10s  gBounds=%s dropped=%s",
					state.name,
					Log.formatBigEngineering(statez.zSumBounds),
					Log.formatBigEngineering(statez.zSumDropped)
				));
			}
		}
		return buf.toString();
	}
}
