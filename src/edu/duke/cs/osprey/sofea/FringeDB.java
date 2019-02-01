package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.ByteBufferInputStream;
import edu.duke.cs.osprey.tools.ByteBufferOutputStream;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.nio.ByteBuffer;
import java.util.Arrays;
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


	private class IOState {

		final BigDecimal[] readZSumMax = new BigDecimal[confSpace.states.size()];
		final BigDecimal[] writeZSumMax = new BigDecimal[confSpace.states.size()];
		long readIndex = 0;
		long numToRead = 0;
		long writeIndex = 0;
		long numWritten = 0;

		IOState() {
			Arrays.fill(readZSumMax, null);
			Arrays.fill(writeZSumMax, null);
		}

		void copyTo(IOState other) {
			System.arraycopy(this.readZSumMax, 0, other.readZSumMax, 0, confSpace.states.size());
			System.arraycopy(this.writeZSumMax, 0, other.writeZSumMax, 0, confSpace.states.size());
			other.readIndex = this.readIndex;
			other.numToRead = this.numToRead;
			other.writeIndex = this.writeIndex;
			other.numWritten = this.numWritten;
		}

		IOState copy() {
			IOState copy = new IOState();
			this.copyTo(copy);
			return copy;
		}

		void advanceRead(int count) {
			readIndex = advanceEntryIndex(readIndex, count);
			numToRead -= count;
		}

		void advanceWrite(int count) {
			writeIndex = advanceEntryIndex(writeIndex, count);
			numWritten += count;
		}

		long advanceEntryIndex(long index, int count) {
			return (index + count) % maxNumEntries;
		}
	}

	public final MultiStateConfSpace confSpace;
	public final File file;

	private final RandomAccessFile io;
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
			iostate.readIndex = io.readLong();
			iostate.numToRead = io.readLong();
			iostate.writeIndex = io.readLong();
			iostate.numWritten = io.readLong();

			// skip to alignment boundary
			for (int i=52; i<64; i++) {
				io.readByte();
			}
			assert (io.getFilePointer() % 32 == 0);

			// read the z stats
			posZStats = io.getFilePointer();
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.readZSumMax[state.index] = bdio.read(io);
			}
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.writeZSumMax[state.index] = bdio.read(io);
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
			maxNumEntries = (file.length() - posEntries)/entryBytes;

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
		return stateEncoding.numBytes + confBytes + bdio.numBytes*4;
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
	public BigDecimal getZSumMax(MultiStateConfSpace.State state) {
		return iostate.readZSumMax[state.index];
	}

	public class Transaction {

		private final IOState iostate = FringeDB.this.iostate.copy();

		private MultiStateConfSpace.State state;
		private int[] conf;
		private BigDecimalBounds zSumBounds;
		private BigDecimalBounds zPathHeadBounds;

		private final ByteBuffer readBuf = ByteBuffer.allocate(1024*1024);
		private final DataInput readIn = new DataInputStream(new ByteBufferInputStream(readBuf));
		private final int maxReadEntries = readBuf.capacity()/entryBytes;

		private final ByteBuffer writeBuf = ByteBuffer.allocate(1024*1024);
		private final DataOutput writeOut = new DataOutputStream(new ByteBufferOutputStream(writeBuf));
		private final int maxWrittenEntries = writeBuf.capacity()/entryBytes;
		private int writtenEntries = 0;

		private Transaction() {
			// keep the constructor private

			readBuf.limit(0);
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

				// if the read buffer is empty, fill it up
				if (readBuf.position() == readBuf.limit()) {
					int numToRead = Math.min(maxReadEntries, (int)Math.min(iostate.numToRead, Integer.MAX_VALUE));
					io.seek(posEntries + iostate.readIndex*entryBytes);
					io.readFully(readBuf.array(), 0, numToRead*entryBytes);
					readBuf.position(0);
					readBuf.limit(numToRead*entryBytes);
				}

				// read the next entry out of the read buffer

				state = confSpace.states.get(stateEncoding.read(readIn));

				// read the conf, and undo the shift
				conf = new int[state.confSpace.positions.size()];
				for (int i=0; i<conf.length; i++) {
					conf[i] = confEncoding.read(readIn) - 1;
				}

				// pad the conf if needed
				int numBytesRead = confEncoding.numBytes*conf.length;
				for (int i=numBytesRead; i<confBytes; i++) {
					readIn.readByte();
				}

				zSumBounds = new BigDecimalBounds(
					bdio.read(readIn),
					bdio.read(readIn)
				);
				zPathHeadBounds = new BigDecimalBounds(
					bdio.read(readIn),
					bdio.read(readIn)
				);

			} catch (IOException ex) {
				throw new RuntimeException("can't advance to next fringe node", ex);
			}

			// advance the read offset, wrapping if needed
			iostate.advanceRead(1);
		}

		public MultiStateConfSpace.State state() {
			return state;
		}

		public int[] conf() {
			return conf;
		}

		public BigDecimalBounds zSumBounds() {
			return zSumBounds;
		}

		public BigDecimalBounds zPathHeadBounds() {
			return zPathHeadBounds;
		}

		private void updateZMax(int stateIndex, BigDecimal val) {
			if (iostate.writeZSumMax[stateIndex] == null || MathTools.isGreaterThan(val, iostate.writeZSumMax[stateIndex])) {
				iostate.writeZSumMax[stateIndex] = val;
			}
		}

		/**
		 * the total number of nodes that can fit in the transaction write buffer
		 */
		public int maxWriteBufferNodes() {
			return maxWrittenEntries;
		}

		/**
		 * Is there enough room in the transaction write buffer to add more nodes?
		 */
		public boolean txHasRoomFor(int count) {
			return writtenEntries + count <= maxWrittenEntries;
		}

		private void writeEntry(int stateIndex, int[] conf, BigDecimalBounds zSumBounds, BigDecimalBounds zPathHeadBounds, DataOutput out) {
			try {

				stateEncoding.write(out, stateIndex);

				// write the conf, but +1 each rc so the min value is 0 instead of -1
				for (int rc : conf) {
					confEncoding.write(out, rc + 1);
				}

				// pad the conf if needed
				int numBytesWritten = confEncoding.numBytes*conf.length;
				assert (numBytesWritten <= confBytes);
				for (int i=numBytesWritten; i<confBytes; i++) {
					out.writeByte(0);
				}

				bdio.write(out, zSumBounds.lower);
				bdio.write(out, zSumBounds.upper);
				bdio.write(out, zPathHeadBounds.lower);
				bdio.write(out, zPathHeadBounds.upper);

			} catch (IOException ex) {
				throw new RuntimeException("can't write fringe node", ex);
			}
		}

		/**
		 * add a node to the transaction write buffer with the same state as the last-read node in the sweep
		 */
		public void writeRootNode(MultiStateConfSpace.State state, BigDecimalBounds zSumBounds, BigDecimalBounds zPathHeadBounds) {

			if (!txHasRoomFor(1)) {
				throw new IllegalStateException("transaction write buffer has no more room for nodes");
			}

			writeEntry(state.index, Conf.make(state.confSpace), zSumBounds, zPathHeadBounds, writeOut);
			writtenEntries++;
			updateZMax(state.index, zSumBounds.upper);
		}

		/**
		 * add a node to the transaction write buffer with the same state as the last-read node in the sweep
		 */
		public void writeReplacementNode(MultiStateConfSpace.State state, int[] conf, BigDecimalBounds zSumBounds, BigDecimalBounds zPathHeadBounds) {

			if (!txHasRoomFor(1)) {
				throw new IllegalStateException("transaction write buffer has no more room for nodes");
			}

			writeEntry(state.index, conf, zSumBounds, zPathHeadBounds, writeOut);
			writtenEntries++;
			updateZMax(state.index, zSumBounds.upper);
		}

		/**
		 * Is there enough room in the database to add more nodes?
		 */
		public boolean dbHasRoomFor(int count) {
			long usedEntries = iostate.numToRead + iostate.numWritten + writtenEntries;
			long freeEntries = maxNumEntries - usedEntries;
			return count <= freeEntries;
		}

		/**
		 * Is there enough room in the database to add the nodes in this transaction?
		 */
		public boolean dbHasRoomForCommit() {
			return dbHasRoomFor(0);
		}

		/**
		 * Flushes all pending writes to the database file.
		 * All writes are flushed to the underlying storage device by the time this method returns.
		 */
		public void commit() {

			// short circuit
			if (writtenEntries <= 0 && iostate.numToRead == FringeDB.this.iostate.numToRead) {
				return;
			}

			if (!dbHasRoomForCommit()) {
				throw new IllegalStateException("transaction too big to commit");
			}

			try {

				// write the replacement entries
				writeBuf.flip();
				io.seek(posEntries + iostate.writeIndex*entryBytes);
				io.write(writeBuf.array(), 0, writtenEntries*entryBytes);
				writeBuf.clear();
				iostate.advanceWrite(writtenEntries);
				writtenEntries = 0;

				// write zSumMax
				io.seek(posZStats + bdio.numBytes*confSpace.states.size());
				for (MultiStateConfSpace.State state : confSpace.states) {
					bdio.write(io, iostate.writeZSumMax[state.index]);
				}

				// persist io state
				io.seek(posIOState);
				io.writeLong(iostate.readIndex);
				io.writeLong(iostate.numToRead);
				io.writeLong(iostate.writeIndex);
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
				iostate.readZSumMax[state.index] = iostate.writeZSumMax[state.index];
				bdio.write(io, iostate.readZSumMax[state.index]);
			}
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.writeZSumMax[state.index] = null;
				bdio.write(io, iostate.writeZSumMax[state.index]);
			}

			// persist read state
			io.seek(posIOState);
			io.writeLong(iostate.readIndex);
			io.writeLong(iostate.numToRead);
			io.writeLong(iostate.writeIndex);
			io.writeLong(iostate.numWritten);

			// flush changes to storage
			io.getChannel().force(false);

		} catch (IOException ex) {
			throw new RuntimeException("finish failed", ex);
		}
	}
}
