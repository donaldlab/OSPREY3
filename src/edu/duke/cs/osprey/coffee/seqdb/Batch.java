package edu.duke.cs.osprey.coffee.seqdb;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Consumer;


public class Batch {

	public final SeqDB seqdb;

	final Map<Sequence,SeqInfo> sequencedSums = new HashMap<>();
	final Map<Integer,StateZ> unsequencedSums = new HashMap<>();

	private boolean isEmpty = true;

	Batch(SeqDB seqdb) {
		this.seqdb = seqdb;
	}

	private void update(MultiStateConfSpace.State state, Sequence seq, Consumer<StateZ> f) {

		if (state.isSequenced) {

			assert (seq != null);

			// get the batch seq info, or empty sums
			SeqInfo seqInfo = sequencedSums.get(seq);
			if (seqInfo == null) {
				seqInfo = SeqInfo.makeZero(seqdb.confSpace);
			}

			f.accept(seqInfo.get(state));
			sequencedSums.put(seq, seqInfo);

		} else {

			assert (seq == null);

			// get the batch sum, or empty
			var statez = unsequencedSums.get(state.unsequencedIndex);
			if (statez == null) {
				statez = StateZ.makeZero();
			}

			f.accept(statez);
			unsequencedSums.put(state.unsequencedIndex, statez);
		}

		isEmpty = false;
	}

	public void addZConf(MultiStateConfSpace.State state, Sequence seq, BigDecimal zConf, BigExp zSumUpper) {

		if (!MathTools.isFinite(zConf) || !zSumUpper.isFinite()) {
			throw new IllegalArgumentException("Z must be finite: " + zConf + ", " + zSumUpper);
		}

		update(state, seq, statez -> {
			statez.zSumBounds.lower = seqdb.bigMath()
				.set(statez.zSumBounds.lower)
				.add(zConf)
				.get();
			statez.zSumBounds.upper = seqdb.bigMath()
				.set(statez.zSumBounds.upper)
				.add(zConf)
				.sub(zSumUpper)
				.get();
		});
	}

	public void addZSumUpper(MultiStateConfSpace.State state, Sequence seq, BigExp zSumUpper) {

		if (!zSumUpper.isFinite()) {
			throw new IllegalArgumentException("Z must be finite: " + zSumUpper);
		}

		update(state, seq, statez ->
			statez.zSumBounds.upper = seqdb.bigMath()
				.set(statez.zSumBounds.upper)
				.add(zSumUpper)
				.get()
		);
	}

	public void subZSumUpper(MultiStateConfSpace.State state, Sequence seq, BigExp zSumUpper) {

		if (!zSumUpper.isFinite()) {
			throw new IllegalArgumentException("Z must be finite: " + zSumUpper);
		}

		update(state, seq, statez ->
			statez.zSumBounds.upper = seqdb.bigMath()
				.set(statez.zSumBounds.upper)
				.sub(zSumUpper)
				.get()
		);
	}

	public void drop(MultiStateConfSpace.State state, Sequence seq, BigExp zSumUpper) {

		if (!zSumUpper.isFinite()) {
			throw new IllegalArgumentException("Z must be finite: " + zSumUpper);
		}

		update(state, seq, statez ->
			statez.zSumDropped = seqdb.bigMath()
				.set(statez.zSumDropped)
				.add(zSumUpper)
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

		if (seqdb.member.isDirector()) {
			// save locally
			seqdb.commitBatch(op);
		} else {
			// relay batch save to the driver member
			seqdb.member.sendTo(op, seqdb.member.directorAddress());
		}
	}

	@Override
	public String toString() {
		var buf = new StringBuilder();
		unsequencedSums.forEach((statei, statez) -> {
			var state = seqdb.confSpace.states.get(statei);
			buf.append(String.format("%10s   bounds %s   dropped %s\n",
				state.name,
				Log.formatBigEngineering(statez.zSumBounds),
				Log.formatBigEngineering(statez.zSumDropped)
			));
		});
		sequencedSums.forEach((seq, seqInfo) -> {
			for (var state : seqdb.confSpace.sequencedStates) {
				var statez = seqInfo.get(state);
				if (!statez.zSumBounds.isZero()) {
					buf.append(String.format("%10s   bounds %s   dropped %s   [%s]\n",
						state.name,
						Log.formatBigEngineering(statez.zSumBounds),
						Log.formatBigEngineering(statez.zSumDropped),
						seq
					));
				}
			}
		});
		return buf.toString();
	}
}
