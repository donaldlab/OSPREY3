package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;


import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.kstar.NewalgLab;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

public class FringeDB {

	private static class Entry {

		final int[] conf;
		final BigDecimalBounds bounds;
		final BigDecimal zpath;

		public Entry(int[] conf, BigDecimalBounds bounds, BigDecimal zpath) {
			this.conf = conf;
			this.bounds = bounds;
			this.zpath = zpath;
		}
	}

	public final SimpleConfSpace confSpace;

	private List<Entry> entriesRead = new ArrayList<>();
	private List<Entry> entriesWrite = new ArrayList<>();
	private BigDecimal zmax = BigDecimal.ZERO;

	public FringeDB(SimpleConfSpace confSpace) {
		this.confSpace = confSpace;
	}

	// TODO: actually persist to disk

	public void add(ConfIndex index, BigDecimalBounds bounds, BigDecimal zpath) {

		// write the entry
		entriesWrite.add(new Entry(Conf.make(index), bounds, zpath));

		// update zmax
		if (MathTools.isGreaterThan(bounds.upper, zmax)) {
			zmax = bounds.upper;
		}
	}

	public BigDecimal getZMax() {
		return zmax;
	}

	public interface Listener {
		void tree(ConfIndex index, BigDecimalBounds bounds, BigDecimal zpath);
	}

	private void swap() {

		// swap the lists
		List<Entry> temp = entriesRead;
		entriesRead = entriesWrite;
		entriesWrite = temp;

		zmax = BigDecimal.ZERO;
	}

	public int writeSize() {
		return entriesWrite.size();
	}

	public void sweep(Listener listener) {

		// assume the things we want to read were just written to the write list
		swap();

		// clear the write list
		entriesWrite.clear();

		// sweep over the read list
		ConfIndex index = new ConfIndex(confSpace.positions.size());
		for (Entry entry : entriesRead) {
			Conf.index(entry.conf, index);
			listener.tree(index, entry.bounds, entry.zpath);
		}
	}

	// TEMP
	public void dump() {
		log("Fringe DB:");
		for (Entry entry : entriesWrite) {
			log("%s  %s", NewalgLab.dump(entry.bounds), Conf.toString(entry.conf));
		}
	}
}
