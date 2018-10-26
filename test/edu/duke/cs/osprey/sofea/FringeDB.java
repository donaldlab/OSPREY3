package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;


import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.kstar.NewalgLab;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class FringeDB {

	private static class Entry {

		final int stateIndex;
		final int[] conf;
		final BigDecimalBounds bounds;
		final BigDecimal zpath;

		public Entry(int stateIndex, int[] conf, BigDecimalBounds bounds, BigDecimal zpath) {
			this.stateIndex = stateIndex;
			this.conf = conf;
			this.bounds = bounds;
			this.zpath = zpath;
		}
	}

	public final MultiStateConfSpace confSpace;

	private List<Entry> entriesRead = new ArrayList<>();
	private List<Entry> entriesWrite = new ArrayList<>();
	private BigDecimal[] zmax;

	public FringeDB(MultiStateConfSpace confSpace) {
		this.confSpace = confSpace;
		this.zmax = new BigDecimal[confSpace.states.size()];
		Arrays.fill(this.zmax, BigDecimal.ZERO);
	}

	// TODO: actually persist to disk

	public void add(MultiStateConfSpace.State state, ConfIndex index, BigDecimalBounds bounds, BigDecimal zpath) {

		// write the entry
		entriesWrite.add(new Entry(state.index, Conf.make(index), bounds, zpath));

		// update zmax
		if (MathTools.isGreaterThan(bounds.upper, zmax[state.index])) {
			zmax[state.index] = bounds.upper;
		}
	}

	public BigDecimal getZMax(MultiStateConfSpace.State state) {
		return zmax[state.index];
	}

	public interface Listener {
		void tree(MultiStateConfSpace.State state, ConfIndex index, BigDecimalBounds bounds, BigDecimal zpath);
	}

	private void swap() {

		// swap the lists
		List<Entry> temp = entriesRead;
		entriesRead = entriesWrite;
		entriesWrite = temp;

		// reset the max Z values
		Arrays.fill(this.zmax, BigDecimal.ZERO);
	}

	public int writeSize() {
		return entriesWrite.size();
	}

	public void sweep(Listener listener) {

		// assume the things we want to read were just written to the write list
		swap();

		// clear the write list
		entriesWrite.clear();

		// allocate a conf index for each state
		List<ConfIndex> indices = confSpace.states.stream()
			.map(state -> new ConfIndex(state.confSpace.positions.size()))
			.collect(Collectors.toList());

		// sweep over the read list
		for (Entry entry : entriesRead) {
			ConfIndex index = indices.get(entry.stateIndex);
			Conf.index(entry.conf, index);
			MultiStateConfSpace.State state = confSpace.states.get(entry.stateIndex);
			listener.tree(state, index, entry.bounds, entry.zpath);
		}
	}

	// TEMP
	public void dump() {
		log("Fringe DB:");
		for (Entry entry : entriesWrite) {
			MultiStateConfSpace.State state = confSpace.states.get(entry.stateIndex);
			log("%s  %s %s", NewalgLab.dump(entry.bounds), state.name, Conf.toString(entry.conf));
		}
	}
}
