package edu.duke.cs.osprey.coffee.seqdb;

import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.coffee.Serializers;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.*;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import org.mapdb.*;

import java.io.File;
import java.math.BigDecimal;
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
	final HTreeMap<Integer,BigDecimalBounds> unsequencedSums;

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
				.valueSerializer(new MapDBTools.BigDecimalBoundsSerializer())
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

		for (var sum : op.sequencedSums) {

			// convert the sum to a SeqInfo
			SeqInfo seqInfo = new SeqInfo(confSpace);
			System.arraycopy(sum.boundsByState, 0, seqInfo.zSumBounds, 0, confSpace.sequencedStates.size());

			// combine with the old sums if needed
			SeqInfo oldSeqInfo = sequencedSums.get(sum.seq);
			if (oldSeqInfo != null) {
				for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
					MathTools.BigDecimalBounds zSum = sum.boundsByState[state.sequencedIndex];
					combineSums(zSum, oldSeqInfo.zSumBounds[state.sequencedIndex]);
					fixRoundoffError(zSum);
				}
			}

			sequencedSums.put(sum.seq, seqInfo);
		}

		for (var sum : op.unsequencedSums) {

			BigDecimalBounds zSum = sum.bounds;

			// combine with the old sum if needed
			MathTools.BigDecimalBounds oldZSum = unsequencedSums.get(sum.unsequencedIndex);
			if (oldZSum != null) {
				combineSums(zSum, oldZSum);
				fixRoundoffError(zSum);
			}

			unsequencedSums.put(sum.unsequencedIndex, zSum);
		}

		db.commit();
	}

	private void combineSums(MathTools.BigDecimalBounds sum, MathTools.BigDecimalBounds oldSum) {
		sum.upper = bigMath()
			.set(sum.upper)
			.add(oldSum.upper)
			.get();
		sum.lower = bigMath()
			.set(sum.lower)
			.add(oldSum.lower)
			.get();
	}

	private void fixRoundoffError(MathTools.BigDecimalBounds z) {
		// trust the lower bound more, since it's based on minimizations
		if (!z.isValid()) {
			// TODO: throw an Exception if the error is bigger than what we'd expect from roundoff?
			z.upper = z.lower;
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
	BigDecimalBounds getSum(MultiStateConfSpace.State state) {

		checkDirector();

		return unsequencedSums.get(state.unsequencedIndex);
	}

	/**
	 * Returns the current Z bounds for the unsequenced state.
	 */
	public BigDecimalBounds boundsUnsequenced(MultiStateConfSpace.State state) {

		checkDirector();

		BigDecimalBounds z = getSum(state);
		if (z == null) {
			z = new BigDecimalBounds(BigDecimal.ZERO, MathTools.BigPositiveInfinity);
		}
		return z;
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
	public SeqInfo boundsUnexplored(Sequence seq) {

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
		seqInfo = new SeqInfo(confSpace);
		seqInfo.setUnknown();
		return seqInfo;
	}

	/**
	 * Returns the current sequenced state Z bounds for the fully-assigned sequence.
	 */
	public SeqInfo bounds(Sequence seq) {

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
		seqInfo = new SeqInfo(confSpace);
		seqInfo.setUnknown();
		return seqInfo;
	}

	/**
	 * returns Z bounds for all sequences
	 * returns bounds for both full and partial sequences
	 */
	public Iterable<Map.Entry<Sequence,SeqInfo>> boundsSequenced() {

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
					seqInfo.zSumBounds[state.sequencedIndex].upper = bigMath()
						.set(seqInfo.zSumBounds[state.sequencedIndex].upper)
						.add(parentSeqInfo.zSumBounds[state.sequencedIndex].upper)
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
			BigDecimalBounds zSumBounds = boundsUnsequenced(state);
			buf.append(String.format("\n%10s  zSumBounds=%s w=%s",
				state.name,
				Log.formatBigLn(zSumBounds),
				Log.formatBigLn(zSumBounds.size(mathContext))
			));
		}
		buf.append("\nSequenced");
		for (Map.Entry<Sequence,SeqInfo> entry : boundsSequenced()) {
			Sequence seq = entry.getKey();
			SeqInfo seqInfo = entry.getValue();
			buf.append(String.format("\n[%s]", seq));
			for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
				buf.append(String.format("\n%10s  zSumBounds=%s w=%s",
					state.name,
					Log.formatBigLn(seqInfo.get(state)),
					Log.formatBigLn(seqInfo.get(state).size(mathContext))
				));
			}
		}
		return buf.toString();
	}
}
