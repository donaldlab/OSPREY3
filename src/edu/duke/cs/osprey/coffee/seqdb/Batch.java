package edu.duke.cs.osprey.coffee.seqdb;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

import java.util.HashMap;
import java.util.Map;
import java.util.function.Consumer;


public class Batch {

	public final SeqDB seqdb;

	final Map<Sequence,SeqInfo> sequencedSums = new HashMap<>();
	final Map<Integer,BigDecimalBounds> unsequencedSums = new HashMap<>();

	private boolean isEmpty = true;

	Batch(SeqDB seqdb) {
		this.seqdb = seqdb;
	}

	private void updateZSumBounds(MultiStateConfSpace.State state, Sequence seq, Consumer<BigDecimalBounds> f) {

		if (state.isSequenced) {

			assert (seq != null);

			// get the batch seq info, or empty sums
			SeqInfo seqInfo = sequencedSums.get(seq);
			if (seqInfo == null) {
				seqInfo = new SeqInfo(seqdb.confSpace);
				seqInfo.setEmpty();
			}

			f.accept(seqInfo.get(state));
			sequencedSums.put(seq, seqInfo);

		} else {

			assert (seq == null);

			// get the batch sum, or empty
			BigDecimalBounds sum = unsequencedSums.get(state.unsequencedIndex);
			if (sum == null) {
				sum = BigDecimalBounds.makeZero();
			}

			f.accept(sum);
			unsequencedSums.put(state.unsequencedIndex, sum);
		}

		isEmpty = false;
	}

	public void addZConf(MultiStateConfSpace.State state, Sequence seq, BigExp zConf, BigExp zSumUpper) {

		if (!zConf.isFinite() || !zSumUpper.isFinite()) {
			throw new IllegalArgumentException("Z must be finite: " + zConf + ", " + zSumUpper);
		}

		updateZSumBounds(state, seq, sum -> {
			sum.lower = seqdb.bigMath()
				.set(sum.lower)
				.add(zConf)
				.get();
			sum.upper = seqdb.bigMath()
				.set(sum.upper)
				.add(zConf)
				.sub(zSumUpper)
				.get();
		});
	}

	public void addZSumUpper(MultiStateConfSpace.State state, Sequence seq, BigExp zSumUpper) {

		if (!zSumUpper.isFinite()) {
			throw new IllegalArgumentException("Z must be finite: " + zSumUpper);
		}

		updateZSumBounds(state, seq, sum ->
			sum.upper = seqdb.bigMath()
				.set(sum.upper)
				.add(zSumUpper)
				.get()
		);
	}

	public void subZSumUpper(MultiStateConfSpace.State state, Sequence seq, BigExp zSumUpper) {

		if (!zSumUpper.isFinite()) {
			throw new IllegalArgumentException("Z must be finite: " + zSumUpper);
		}

		updateZSumBounds(state, seq, sum ->
			sum.upper = seqdb.bigMath()
				.set(sum.upper)
				.sub(zSumUpper)
				.get()
		);
	}

	public boolean isEmpty() {
		return isEmpty;
	}

	public void save() {

		// short circuit
		if (isEmpty) {
			return;
		}

		var op = new SaveOperation(seqdb.confSpace, this);

		if (seqdb.member.isDriver()) {
			// save locally
			seqdb.commitBatch(op);
		} else {
			// relay batch save to the driver member
			seqdb.member.sendTo(op, seqdb.member.driverAddress());
		}
	}
}
