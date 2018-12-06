package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
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

	// TODO: the initial implementation is super slow!!
	// need to optimize... maybe try batching more ops into Sweep.commit()/rollback()?

	static final byte[] Magic = { 'f', 'r', 'i', 'n', 'g', 'e', 'd', 'b' };


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
	public final File file;

	private final RandomAccessFile io;
	private final long iosize;

	private final IntEncoding stateEncoding;
	private final IntEncoding confEncoding;
	private final BigDecimalIO.Fixed bdio;
	private final int confBytes;
	private final int entryBytes;

	private final long posReadWriteState;
	private final long posZStats;
	private final long posEntries;
	private final long maxNumEntries;

	private final BigDecimal[] readZmax;
	private final BigDecimal[] writeZmax;
	private long readOffset = 0;
	private long numToRead = 0;
	private long writeOffset = 0;
	private long numWritten = 0;

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
			posReadWriteState = io.getFilePointer();
			readOffset = io.readLong();
			numToRead = io.readLong();
			writeOffset = io.readLong();
			numWritten = io.readLong();

			// skip to alignment boundary
			for (int i=52; i<64; i++) {
				io.readByte();
			}

			// read the z stats
			posZStats = io.getFilePointer();
			readZmax = new BigDecimal[confSpace.states.size()];
			for (MultiStateConfSpace.State state : confSpace.states) {
				readZmax[state.index] = bdio.read(io);
			}
			writeZmax = new BigDecimal[confSpace.states.size()];
			for (MultiStateConfSpace.State state : confSpace.states) {
				writeZmax[state.index] = bdio.read(io);
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

			bdio.write(io, entry.bounds.lower);
			bdio.write(io, entry.bounds.upper);
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
		return numToRead + numWritten;
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
		return readZmax[state.index];
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
		log("\tpos read/write = %d", posReadWriteState);
		log("\tpos z stats = %d", posZStats);
		log("\tpos entries = %d", posEntries);
		log("\tmax num entries = %d", maxNumEntries);
		log("\tread offset = %d", readOffset);
		log("\tnum to read = %d", numToRead);
		log("\twrite offset = %d", writeOffset);
		log("\tnum written = %d", numWritten);
		for (MultiStateConfSpace.State state : confSpace.states) {
			log("\t\tread Z[%d] = %s", state.index, readZmax[state.index]);
		}
		for (MultiStateConfSpace.State state : confSpace.states) {
			log("\t\twrite Z[%d] = %s", state.index, writeZmax[state.index]);
		}
	}


	public class Roots {

		private final Entry[] entries = new Entry[confSpace.states.size()];

		private Roots() {
			// keep constructor private
		}

		public void set(MultiStateConfSpace.State state, BigDecimalBounds bounds, BigDecimal zpath) {
			entries[state.index] = new Entry(
				state.index,
				Conf.make(state.confSpace),
				bounds,
				zpath
			);
		}

		public boolean hasSpaceToCommit() {
			return entries.length <= maxNumEntries;
		}

		public void commit() {

			if (!hasSpaceToCommit()) {
				throw new IllegalStateException("not enough space to commit");
			}

			try {

				// write the entries
				io.seek(posEntries);
				for (Entry entry : entries) {
					writeEntry(entry);
				}

				// update zmax
				io.seek(posZStats);
				for (int i=0; i<entries.length; i++) {
					readZmax[i] = entries[i].bounds.upper;
					bdio.write(io, readZmax[i]);
				}
				for (int i=0; i<entries.length; i++) {
					writeZmax[i] = null;
					bdio.write(io, writeZmax[i]);
				}

				// update read/write state to prep for next sweep
				readOffset = 0;
				numToRead = entries.length;
				writeOffset = entryBytes*entries.length;
				numWritten = 0;

				// persist read/write state
				io.seek(posReadWriteState);
				io.writeLong(readOffset);
				io.writeLong(numToRead);
				io.writeLong(writeOffset);
				io.writeLong(numWritten);

				// flush changes to storage
				io.getChannel().force(false);

			} catch (IOException ex) {
				throw new RuntimeException("commit failed", ex);
			}
		}
	}

	public Roots roots() {

		if (numToRead > 0) {
			throw new IllegalStateException("already have roots, can't add more");
		}

		return new Roots();
	}


	public class Sweep {

		private final long numNodes = numToRead;
		private final ConfIndex[] indices;
		private BigDecimal zmax;
		private Entry entry;
		private final List<Entry> entries = new ArrayList<>();

		private Sweep() {

			// allocate a conf index for each state as scratch space
			indices = new ConfIndex[confSpace.states.size()];
			for (MultiStateConfSpace.State state : confSpace.states) {
				indices[state.index] = new ConfIndex(state.confSpace.positions.size());
			}
		}

		/** total number of nodes in this sweep */
		public long numNodes() {
			return numNodes;
		}

		/** number of unread nodes left in this sweep */
		public long numNodesRemaining() {
			return numToRead;
		}

		/** returns true if there are no remaining nodes to read, false otherwise */
		public boolean isEmpty() {
			return numToRead <= 0;
		}

		/** reads the node at the current read position */
		public void read() {

			if (numToRead <= 0) {
				throw new NoSuchElementException("out of fringe nodes to read");
			}

			try {
				io.seek(posEntries + readOffset);
				entry = readEntry();
			} catch (IOException ex) {
				throw new RuntimeException("can't advance to next fringe node", ex);
			}

			Conf.index(entry.conf, index());
			zmax = writeZmax[entry.stateIndex];
		}

		public MultiStateConfSpace.State state() {
			return confSpace.states.get(entry.stateIndex);
		}

		public ConfIndex index() {
			return indices[entry.stateIndex];
		}

		public BigDecimalBounds bounds() {
			return entry.bounds;
		}

		public BigDecimal zpath() {
			return entry.zpath;
		}

		/**
		 * add a node to the pending write set
		 * save changes by calling commitAndAdvance()
		 * or discard changes by calling rollbackAndAdvance()
		 */
		public void addChild(ConfIndex index, BigDecimalBounds bounds, BigDecimal zpath) {

			entries.add(new Entry(entry.stateIndex, Conf.make(index), bounds, zpath));

			if (zmax == null || MathTools.isGreaterThan(bounds.upper, zmax)) {
				zmax = bounds.upper;
			}
		}

		public boolean hasSpaceToCommit() {
			long usedEntries = numToRead - 1 + numWritten;
			long freeEntries = maxNumEntries - usedEntries;
			return entries.size() <= freeEntries;
		}

		/**
		 * Removes the read node from the queue and adds the new child nodes to the end of the queue.
		 * All writes are flushed to the underlying storage device by the time this method returns.
		 */
		public void commitAndAdvance() {

			if (!hasSpaceToCommit()) {
				throw new IllegalStateException("not enough space to commit. Must rollback instead");
			}

			try {

				// write the entries
				for (Entry entry : entries) {

					io.seek(posEntries + writeOffset);
					writeEntry(entry);

					// advance the write offset, wrapping if needed
					writeOffset += entryBytes;
					assert (posEntries + writeOffset <= iosize);
					if (posEntries + writeOffset == iosize) {
						writeOffset = 0;
					}
					numWritten++;
				}
				entries.clear();

				// update zmax and persist if needed
				boolean betterZmax = zmax != null
					&& (writeZmax[entry.stateIndex] == null || MathTools.isGreaterThan(zmax, writeZmax[entry.stateIndex]));
				if (betterZmax) {

					writeZmax[entry.stateIndex] = zmax;

					io.seek(posZStats + bdio.numBytes*(confSpace.states.size() + entry.stateIndex));
					bdio.write(io, zmax);
				}

				// advance the read offset, wrapping if needed
				readOffset += entryBytes;
				assert (posEntries + readOffset <= iosize);
				if (posEntries + readOffset == iosize) {
					readOffset = 0;
				}
				numToRead--;

				// persist read/write state
				io.seek(posReadWriteState);
				io.writeLong(readOffset);
				io.writeLong(numToRead);
				io.writeLong(writeOffset);
				io.writeLong(numWritten);

				// flush changes to storage
				io.getChannel().force(false);

			} catch (IOException ex) {
				throw new RuntimeException("commit failed", ex);
			}
		}

		/**
		 * Removes the read node from the queue and adds it to the end of the queue.
		 * Any child nodes are discarded.
		 * All writes are flushed to the underlying storage device by the time this method returns.
		 */
		public void rollbackAndAdvance() {
			try {

				entries.clear();

				if (writeOffset != readOffset) {

					// write the entry back onto the queue
					io.seek(posEntries + writeOffset);
					writeEntry(entry);
				}

				// update zmax and persist if needed
				boolean betterZmax = writeZmax[entry.stateIndex] == null
					|| MathTools.isGreaterThan(entry.bounds.upper, writeZmax[entry.stateIndex]);
				if (betterZmax) {

					writeZmax[entry.stateIndex] = entry.bounds.upper;

					io.seek(posZStats + bdio.numBytes*(confSpace.states.size() + entry.stateIndex));
					bdio.write(io, writeZmax[entry.stateIndex]);
				}

				// advance the read offset, wrapping if needed
				readOffset += entryBytes;
				assert (posEntries + readOffset <= iosize);
				if (posEntries + readOffset == iosize) {
					readOffset = 0;
				}
				numToRead--;

				// advance the write offset, wrapping if needed
				writeOffset += entryBytes;
				assert (posEntries + writeOffset <= iosize);
				if (posEntries + writeOffset == iosize) {
					writeOffset = 0;
				}
				numWritten++;

				// persist read state
				io.seek(posReadWriteState);
				io.writeLong(readOffset);
				io.writeLong(numToRead);
				io.writeLong(writeOffset);
				io.writeLong(numWritten);

				// flush changes to storage
				io.getChannel().force(false);

			} catch (IOException ex) {
				throw new RuntimeException("commit failed", ex);
			}
		}

		/**
		 * Finishes this sweep and prepares for the next sweep.
		 * All writes are flushed to the underlying storage device by the time this method returns.
		 */
		public void finish() {

			if (numToRead > 0) {
				throw new IllegalStateException("sweep not finished, " + numToRead + " nodes left to read");
			}

			try {

				// update the entry counts
				numToRead = numWritten;
				numWritten = 0;

				// update the z stats
				io.seek(posZStats);
				for (MultiStateConfSpace.State state : confSpace.states) {
					readZmax[state.index] = writeZmax[state.index];
					bdio.write(io, readZmax[state.index]);
				}
				for (MultiStateConfSpace.State state : confSpace.states) {
					writeZmax[state.index] = null;
					bdio.write(io, writeZmax[state.index]);
				}

				// persist read state
				io.seek(posReadWriteState);
				io.writeLong(readOffset);
				io.writeLong(numToRead);
				io.writeLong(writeOffset);
				io.writeLong(numWritten);

				// flush changes to storage
				io.getChannel().force(false);

			} catch (IOException ex) {
				throw new RuntimeException("commit failed", ex);
			}
		}
	}

	/** starts a new sweep, or resumes the current sweep */
	public Sweep sweep() {
		return new Sweep();
	}
}
