package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;


/**
 * An efficient and resilient database for storing fringe nodes in bounded space.
 *
 * Uses persistent storage (eg HDD, SSD) to ensure nodes are saved after the program exits,
 * with size limits to prevent the database from consuming all available storage space.
 *
 * Nodes are stored as if in a FIFO queue and are accessed via "sweep" operations:
 *   First, the root nodes are written to the database via the root() method.
 *   Then, the entire queue is read in the same order as written via the sweep() method.
 *   During the sweep, nodes for the next sweep are written if space allows.
 *
 * The FIFO queue is implemented by a fixed-size circular buffer.
 * If the queue fills up, no new child nodes can be written, but
 * the current queue can continue to be swept as many times as needed.
 */
public class FringeDB implements AutoCloseable {

	static final byte[] Magic = { 'f', 'r', 'i', 'n', 'g', 'e', 'd', 'b' };


	private static class Entry {

		final int stateIndex;
		final int[] conf;
		final BigDecimalBounds zbounds;
		final BigDecimal zpath;

		public Entry(int stateIndex, int[] conf, BigDecimalBounds zbounds, BigDecimal zpath) {
			this.stateIndex = stateIndex;
			this.conf = conf;
			this.zbounds = zbounds;
			this.zpath = zpath;
		}
	}

	private class IOState {

		final BigDecimal[] readZmax = new BigDecimal[confSpace.states.size()];
		final BigDecimal[] writeZmax = new BigDecimal[confSpace.states.size()];
		long readOffset = 0;
		long numToRead = 0;
		long writeOffset = 0;
		long numWritten = 0;

		IOState() {
			Arrays.fill(readZmax, null);
			Arrays.fill(writeZmax, null);
		}

		void copyTo(IOState other) {
			System.arraycopy(this.readZmax, 0, other.readZmax, 0, confSpace.states.size());
			System.arraycopy(this.writeZmax, 0, other.writeZmax, 0, confSpace.states.size());
			other.readOffset = this.readOffset;
			other.numToRead = this.numToRead;
			other.writeOffset = this.writeOffset;
			other.numWritten = this.numWritten;
		}

		IOState copy() {
			IOState copy = new IOState();
			this.copyTo(copy);
			return copy;
		}

		void advanceRead() {
			readOffset = advanceEntryOffset(readOffset);
			numToRead--;
		}

		void advanceWrite() {
			writeOffset = advanceEntryOffset(writeOffset);
			numWritten++;
		}

		long advanceEntryOffset(long offset) {
			assert (posEntries + offset <= iosize);
			offset += entryBytes;
			if (posEntries + offset >= iosize) {
				offset = 0;
			}
			return offset;
		}
	}

	public final MultiStateConfSpace confSpace;
	public final File file;

	private final RandomAccessFile io;
	private final long iosize;
	private final IOState iostate;

	private final IntEncoding stateEncoding;
	private final IntEncoding confEncoding;
	private final BigDecimalIO.Fixed bdio;
	private final int confBytes;
	private final int entryBytes;

	private final long posIOState;
	private final long posZStats;
	private final long posEntries;
	private final long maxNumEntries;

	/** create a new fringe node database, reserving the desired spase on the filesystem */
	public static FringeDB create(MultiStateConfSpace confSpace, File file, long sizeBytes, MathContext mathContext) {

		// write the header to a new file
		try (RandomAccessFile io = new RandomAccessFile(file, "rw")) {

			// reserve the requested file size in the file system
			io.setLength(sizeBytes);

			// write the magic number
			for (int i=0; i<8; i++) {
				io.writeByte(Magic[i]);
			}

			// write the version
			io.writeInt(1);

			// how many bytes per big decimal?
			BigDecimalIO.Fixed bdio = new BigDecimalIO.Fixed(mathContext);
			io.writeInt(bdio.numBytes);

			// how many bytes per conf?
			IntEncoding confEncoding = getConfEncoding(confSpace);
			int maxConfSize = confSpace.states.stream()
				.mapToInt(state -> state.confSpace.positions.size())
				.max()
				.orElse(0);
			io.writeInt(confEncoding.numBytes*maxConfSize);

			// write the entry read/write state
			io.writeLong(0);
			io.writeLong(0);
			io.writeLong(0);
			io.writeLong(0);

			// header bytes written: 8 + 4*3 + 4*8 = 52

			// pad to 32 bytes so the rest of the file can be nicely aligned
			for (int i=52; i<64; i++) {
				io.writeByte(0);
			}
			assert (io.getFilePointer() % 32 == 0);

			// init the z stats
			for (MultiStateConfSpace.State state : confSpace.states) {
				bdio.write(io, null);
			}
			for (MultiStateConfSpace.State state : confSpace.states) {
				bdio.write(io, null);
			}

			// pad to 32 bytes
			int numBytesWritten = confSpace.states.size()*bdio.numBytes*2;
			int numBytesNeeded = MathTools.roundUpToMultiple(numBytesWritten, 32);
			for (int i=numBytesWritten; i<numBytesNeeded; i++) {
				io.writeByte(0);
			}
			assert (io.getFilePointer() % 32 == 0);

		} catch (IOException ex) {
			throw new RuntimeException(
				String.format("can't create new FringeDB of size %s (%d B) at %s",
					MathTools.formatBytes(sizeBytes),
					sizeBytes,
					file.getAbsolutePath()
				),
				ex
			);
		}

		return open(confSpace, file);
	}

	private static IntEncoding getConfEncoding(MultiStateConfSpace confSpace) {

		// we're going to do arithmatic with the unassigned value, so make sure it's always -1
		assert (Conf.Unassigned == -1);

		return IntEncoding.get(1 + confSpace.states.stream()
			.mapToInt(state -> state.confSpace.positions.stream()
				.mapToInt(pos -> pos.resConfs.stream()
					.mapToInt(rc -> rc.index)
					.max()
					.orElse(0)
				)
				.max()
				.orElse(0)
			)
			.max()
			.orElse(0)
		);
	}

	/** open an existing fringe node database */
	public static FringeDB open(MultiStateConfSpace confSpace, File file) {
		return new FringeDB(confSpace, file);
	}

	private FringeDB(MultiStateConfSpace confSpace, File file) {

		this.confSpace = confSpace;
		this.file = file;

		iosize = file.length();

		// figure out the encodings
		stateEncoding = IntEncoding.get(confSpace.states.stream()
			.mapToInt(state -> state.index)
			.max()
			.orElse(0)
		);
		confEncoding = getConfEncoding(confSpace);

		// open the file and read the header
		try {

			io = new RandomAccessFile(file, "rw");

			// check the magic number
			for (int i=0; i<8; i++) {
				if (io.readByte() != Magic[i]) {
					throw new IOException("not a fringe db file");
				}
			}

			// check the version
			int version = io.readInt();
			if (version != 1) {
				throw new IOException("unrecognized fringe db version: " + version);
			}

			// read the sizes
			bdio = new BigDecimalIO.Fixed(io.readInt());
			confBytes = io.readInt();
			entryBytes = calcEntrySize();

			// read the read/write state
			posIOState = io.getFilePointer();
			iostate = new IOState();
			iostate.readOffset = io.readLong();
			iostate.numToRead = io.readLong();
			iostate.writeOffset = io.readLong();
			iostate.numWritten = io.readLong();

			// skip to alignment boundary
			for (int i=52; i<64; i++) {
				io.readByte();
			}

			// read the z stats
			posZStats = io.getFilePointer();
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.readZmax[state.index] = bdio.read(io);
			}
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.writeZmax[state.index] = bdio.read(io);
			}

			// pad to 32 bytes
			int numBytesRead = confSpace.states.size()*bdio.numBytes*2;
			int numBytesNeeded = MathTools.roundUpToMultiple(numBytesRead, 32);
			for (int i=numBytesRead; i<numBytesNeeded; i++) {
				io.readByte();
			}
			assert (io.getFilePointer() % 32 == 0);

			posEntries = io.getFilePointer();

			// how many entries can we have?
			maxNumEntries = (iosize - posEntries)/entryBytes;

		} catch (IOException ex) {
			throw new RuntimeException("can't open db file: " + file.getAbsolutePath(), ex);
		}
	}

	@Override
	public void close() {
		try {
			io.close();
		} catch (IOException ex) {
			// don't care
		}
	}

	private int calcEntrySize() {
		return stateEncoding.numBytes + confBytes + bdio.numBytes*3;
	}

	private boolean isAlignedOnEntry()
	throws IOException {
		return (io.getFilePointer() - posEntries) % entryBytes == 0;
	}

	private void writeEntry(Entry entry) {
		try {

			assert (isAlignedOnEntry());

			stateEncoding.write(io, entry.stateIndex);

			// write the conf, but +1 each rc so the min value is 0 instead of -1
			for (int i=0; i<entry.conf.length; i++) {
				confEncoding.write(io, entry.conf[i] + 1);
			}

			// pad the conf if needed
			int numBytesWritten = confEncoding.numBytes*entry.conf.length;
			assert (numBytesWritten <= confBytes);
			for (int i=numBytesWritten; i<confBytes; i++) {
				io.writeByte(0);
			}

			bdio.write(io, entry.zbounds.lower);
			bdio.write(io, entry.zbounds.upper);
			bdio.write(io, entry.zpath);

			assert (isAlignedOnEntry());

		} catch (IOException ex) {
			throw new RuntimeException("can't write fringe node", ex);
		}
	}

	private Entry readEntry() {
		try {

			assert (isAlignedOnEntry());

			int stateIndex = stateEncoding.read(io);

			// read the conf, and undo the shift
			int[] conf = new int[confSpace.states.get(stateIndex).confSpace.positions.size()];
			for (int i=0; i<conf.length; i++) {
				conf[i] = confEncoding.read(io) - 1;
			}

			// pad the conf if needed
			int numBytesRead = confEncoding.numBytes*conf.length;
			for (int i=numBytesRead; i<confBytes; i++) {
				io.readByte();
			}

			BigDecimalBounds bounds = new BigDecimalBounds(
				bdio.read(io),
				bdio.read(io)
			);

			BigDecimal zpath = bdio.read(io);

			assert (isAlignedOnEntry());

			return new Entry(stateIndex, conf, bounds, zpath);

		} catch (IOException ex) {
			throw new RuntimeException("can't read fringe node", ex);
		}
	}

	/**
	 * Returns the total number of nodes in the database
	 */
	public long getNumNodes() {
		return iostate.numToRead + iostate.numWritten;
	}

	/**
	 * Returns the maximum number of fringe nodes that can fit in the database.
	 */
	public long getCapacity() {
		return maxNumEntries;
	}

	/**
	 * Returns true if there are no nodes in the database, false otherwise.
	 */
	public boolean isEmpty() {
		return getNumNodes() <= 0;
	}

	/**
	 * Returns the largest Z value for the nodes to read.
	 * Ignores pending Z values in the written nodes.
	 */
	public BigDecimal getZMax(MultiStateConfSpace.State state) {
		return iostate.readZmax[state.index];
	}

	// TEMP
	public void dump() {
		log("FringeDB:");
		log("\tiosize = %d", iosize);
		log("\tstateEncoding = %s", stateEncoding);
		log("\tconfEncoding = %s", confEncoding);
		log("\tbdio bytes = %d", bdio.numBytes);
		log("\tconf bytes = %d", confBytes);
		log("\tentry bytes = %d", entryBytes);
		log("\tpos read/write = %d", posIOState);
		log("\tpos z stats = %d", posZStats);
		log("\tpos entries = %d", posEntries);
		log("\tmax num entries = %d", maxNumEntries);
		log("\tread offset = %d", iostate.readOffset);
		log("\tnum to read = %d", iostate.numToRead);
		log("\twrite offset = %d", iostate.writeOffset);
		log("\tnum written = %d", iostate.numWritten);
		for (MultiStateConfSpace.State state : confSpace.states) {
			log("\t\tread Z[%d] = %s", state.index, iostate.readZmax[state.index]);
		}
		for (MultiStateConfSpace.State state : confSpace.states) {
			log("\t\twrite Z[%d] = %s", state.index, iostate.writeZmax[state.index]);
		}
	}

	public class Transaction {

		private final IOState iostate = FringeDB.this.iostate.copy();

		private Entry entry;
		private final List<Entry> entriesToWrite = new ArrayList<>();

		private Transaction() {
			// keep the constructor private
		}

		/** number of unread nodes left in this sweep */
		public long numNodesToRead() {
			return iostate.numToRead;
		}

		/** returns true if there are no remaining nodes to read in this sweep, false otherwise */
		public boolean hasNodesToRead() {
			return iostate.numToRead > 0;
		}

		/** reads and removes the node at the head of the queue */
		public void readNode() {

			if (iostate.numToRead <= 0) {
				throw new NoSuchElementException("out of fringe nodes to read");
			}

			try {
				io.seek(posEntries + iostate.readOffset);
				entry = readEntry();
			} catch (IOException ex) {
				throw new RuntimeException("can't advance to next fringe node", ex);
			}

			// advance the read offset, wrapping if needed
			iostate.advanceRead();
		}

		public MultiStateConfSpace.State state() {
			return confSpace.states.get(entry.stateIndex);
		}

		public int[] conf() {
			return entry.conf;
		}

		public BigDecimalBounds zbounds() {
			return entry.zbounds;
		}

		public BigDecimal zpath() {
			return entry.zpath;
		}

		private void updateZMax(int stateIndex, BigDecimal val) {
			if (iostate.writeZmax[stateIndex] == null || MathTools.isGreaterThan(val, iostate.writeZmax[stateIndex])) {
				iostate.writeZmax[stateIndex] = val;
			}
		}

		/**
		 * add a node to the pending write set with the same state as the last-read node in the sweep
		 */
		public void addRootNode(MultiStateConfSpace.State state, BigDecimalBounds zbounds, BigDecimal zpath) {
			int[] conf = Conf.make(state.confSpace);
			entriesToWrite.add(new Entry(state.index, conf, zbounds, zpath));
			updateZMax(state.index, zbounds.upper);
		}

		/**
		 * add a node to the pending write set with the same state as the last-read node in the sweep
		 */
		public void addReplacementNode(int[] conf, BigDecimalBounds zbounds, BigDecimal zpath) {
			entriesToWrite.add(new Entry(entry.stateIndex, conf, zbounds, zpath));
			updateZMax(entry.stateIndex, zbounds.upper);
		}

		/**
		 * Is there enough room in the database to add more nodes?
		 */
		public boolean hasRoomFor(int count) {
			long usedEntries = iostate.numToRead + iostate.numWritten + entriesToWrite.size();
			long freeEntries = maxNumEntries - usedEntries;
			return count <= freeEntries;
		}

		/**
		 * Is there enough room in the database to add the nodes in this transaction?
		 */
		public boolean hasRoomForCommit() {
			return hasRoomFor(0);
		}

		/**
		 * Flushes all pending writes to the database file.
		 * All writes are flushed to the underlying storage device by the time this method returns.
		 */
		public void commit() {

			if (!hasRoomForCommit()) {
				throw new IllegalStateException("transaction too big to commit");
			}

			try {

				// write the replacement entries
				for (Entry entry : entriesToWrite) {

					io.seek(posEntries + iostate.writeOffset);
					writeEntry(entry);
					iostate.advanceWrite();
				}
				entriesToWrite.clear();

				// write zmax
				io.seek(posZStats + bdio.numBytes*confSpace.states.size());
				for (MultiStateConfSpace.State state : confSpace.states) {
					bdio.write(io, iostate.writeZmax[state.index]);
				}

				// persist io state
				io.seek(posIOState);
				io.writeLong(iostate.readOffset);
				io.writeLong(iostate.numToRead);
				io.writeLong(iostate.writeOffset);
				io.writeLong(iostate.numWritten);

				// copy io state outside of transaction
				iostate.copyTo(FringeDB.this.iostate);

				// flush changes to storage
				io.getChannel().force(false);

			} catch (IOException ex) {
				throw new RuntimeException("commit failed", ex);
			}
		}
	}

	/** starts a new transaction in the current sweep */
	public Transaction transaction() {
		return new Transaction();
	}

	/**
	 * Finishes this sweep and prepares for the next sweep.
	 * All writes are flushed to the underlying storage device by the time this method returns.
	 */
	public void finishSweep() {

		if (iostate.numToRead > 0) {
			throw new IllegalStateException("sweep not finished, " + iostate.numToRead + " nodes left to read");
		}

		try {

			// update the entry counts
			iostate.numToRead = iostate.numWritten;
			iostate.numWritten = 0;

			// update the z stats
			io.seek(posZStats);
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.readZmax[state.index] = iostate.writeZmax[state.index];
				bdio.write(io, iostate.readZmax[state.index]);
			}
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.writeZmax[state.index] = null;
				bdio.write(io, iostate.writeZmax[state.index]);
			}

			// persist read state
			io.seek(posIOState);
			io.writeLong(iostate.readOffset);
			io.writeLong(iostate.numToRead);
			io.writeLong(iostate.writeOffset);
			io.writeLong(iostate.numWritten);

			// flush changes to storage
			io.getChannel().force(false);

		} catch (IOException ex) {
			throw new RuntimeException("finish failed", ex);
		}
	}
}
