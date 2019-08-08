/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.*;

import java.io.*;
import java.nio.ByteBuffer;
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

		final BigExp[] readZSumMax = new BigExp[confSpace.states.size()];
		final BigExp[] writeZSumMax = new BigExp[confSpace.states.size()];
		long readIndex = 0;
		long numToRead = 0;
		long writeIndex = 0;
		long numWritten = 0;

		IOState() {
			for (MultiStateConfSpace.State state : confSpace.states) {
				readZSumMax[state.index] = new BigExp(Double.NaN);
				writeZSumMax[state.index] = new BigExp(Double.NaN);
			}
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
	private final int confBytes;
	private final int entryBytes;

	private final long posIOState;
	private final long posZStats;
	private final long posEntries;
	private final long maxNumEntries;

	/** create a new fringe node database, reserving the desired spase on the filesystem */
	public static FringeDB create(MultiStateConfSpace confSpace, File file, long sizeBytes) {

		if (sizeBytes <= 0) {
			throw new IllegalArgumentException("invalid FringeDB size: " + sizeBytes + " bytes");
		}

		// write the header to a new file
		try (RandomAccessFile io = new RandomAccessFile(file, "rw")) {

			// reserve the requested file size in the file system
			io.setLength(sizeBytes);

			// write the magic number
			for (int i=0; i<8; i++) {
				io.writeByte(Magic[i]);
			}

			// write the version
			io.writeInt(2);

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

			// header bytes written: 8 + 2*4 + 4*8 = 48

			// pad to 32 bytes so the rest of the file can be nicely aligned
			for (int i=48; i<64; i++) {
				io.writeByte(0);
			}
			assert (io.getFilePointer() % 32 == 0);

			// init the z stats
			for (int i=0; i<confSpace.states.size()*2; i++) {
				new BigExp(Double.NaN).writeTo(io);
			}

			// pad to 32 bytes
			int numBytesWritten = confSpace.states.size()*BigExp.NumBytes*2;
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
			if (version != 2) {
				throw new IOException("unrecognized fringe db version: " + version);
			}

			// read the sizes
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
			for (int i=48; i<64; i++) {
				io.readByte();
			}
			assert (io.getFilePointer() % 32 == 0);

			// read the z stats
			posZStats = io.getFilePointer();
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.readZSumMax[state.index].readFrom(io);
			}
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.writeZSumMax[state.index].readFrom(io);
			}

			// pad to 32 bytes
			int numBytesRead = confSpace.states.size()*BigExp.NumBytes*2;
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
		return stateEncoding.numBytes + confBytes + BigExp.NumBytes;
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

	/** number of unread nodes left in this sweep */
	public long numNodesToRead() {
		return iostate.numToRead;
	}

	/** returns true if there are no remaining nodes to read in this sweep, false otherwise */
	public boolean hasNodesToRead() {
		return iostate.numToRead > 0;
	}

	/**
	 * Returns the largest Z value for the nodes to read.
	 * Ignores pending Z values in the written nodes.
	 */
	public BigExp getZSumMax(MultiStateConfSpace.State state) {
		return iostate.readZSumMax[state.index];
	}

	public class Transaction {

		private final IOState iostate = FringeDB.this.iostate.copy();

		private MultiStateConfSpace.State state;
		private int[] conf;
		private BigExp zSumUpper;

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

					// what's the most entries we could read without wrapping around?
					long numToRead = maxNumEntries - iostate.readIndex;

					// but don't read more entries than are waiting
					numToRead = Math.min(numToRead, iostate.numToRead);

					// but don't read more than the buffer
					numToRead = Math.min(numToRead, maxReadEntries);

					try {

						int readSize = (int)(numToRead*entryBytes);
						io.seek(posEntries + iostate.readIndex*entryBytes);
						io.readFully(readBuf.array(), 0, readSize);
						readBuf.position(0);
						readBuf.limit(readSize);

					} catch (IOException ex) {

						// this actually happened, so add more debugging info to the exception
						throw new IOException(
							String.format(
								"can't fill read buffer."
									+ "\n\tmaxReadEntries=%d"
									+ "\n\tiostate.numToRead=%d"
									+ "\n\tmaxNumToRead=%d"
									+ "\n\tnumToRead=%d"
									+ "\n\tposEntries=%d"
									+ "\n\tiostate.readIndex=%d"
									+ "\n\tentryBytes=%d"
									+ "\n\tseek pos=%d"
									+ "\n\tread len=%d"
									+ "\n\tend read pos=%d"
									+ "\n\tcapacity=%d",
								maxReadEntries,
								iostate.numToRead,
								maxNumEntries - iostate.readIndex,
								numToRead,
								posEntries,
								iostate.readIndex,
								entryBytes,
								posEntries + iostate.readIndex*entryBytes,
								numToRead*entryBytes,
								posEntries + iostate.readIndex*entryBytes + numToRead*entryBytes,
								io.length()
							),
							ex
						);
					}
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

				zSumUpper = BigExp.read(readIn);

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

		public BigExp zSumUpper() {
			return zSumUpper;
		}

		private void updateZMax(int stateIndex, BigExp val) {
			if (iostate.writeZSumMax[stateIndex].isNaN() || val.greaterThan(iostate.writeZSumMax[stateIndex])) {
				iostate.writeZSumMax[stateIndex].set(val);
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

		private void writeEntry(int stateIndex, int[] conf, BigExp zSumUpper, DataOutput out) {
			try {

				// write the state index
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

				zSumUpper.writeTo(out);

			} catch (IOException ex) {
				throw new RuntimeException("can't write fringe node", ex);
			}
		}

		/**
		 * add a node to the transaction write buffer with the same state as the last-read node in the sweep
		 */
		public void writeRootNode(MultiStateConfSpace.State state, BigExp zSumUpper) {

			if (!txHasRoomFor(1)) {
				throw new IllegalStateException("transaction write buffer has no more room for nodes");
			}

			writeEntry(state.index, Conf.make(state.confSpace), zSumUpper, writeOut);
			writtenEntries++;
			updateZMax(state.index, zSumUpper);
		}

		/**
		 * add a node to the transaction write buffer with the same state as the last-read node in the sweep
		 */
		public void writeReplacementNode(MultiStateConfSpace.State state, int[] conf, BigExp zSumUpper) {

			if (!txHasRoomFor(1)) {
				throw new IllegalStateException("transaction write buffer has no more room for nodes");
			}

			writeEntry(state.index, conf, zSumUpper, writeOut);
			writtenEntries++;
			updateZMax(state.index, zSumUpper);
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
				long maxNumToWriteAtOnce = maxNumEntries - iostate.writeIndex;
				if (writtenEntries <= maxNumToWriteAtOnce) {

					// write it all in one pass
					io.seek(posEntries + iostate.writeIndex*entryBytes);
					io.write(writeBuf.array(), 0, writtenEntries*entryBytes);
					iostate.advanceWrite(writtenEntries);

				} else {

					// write in two passes
					int numWrittenPass1 = (int)maxNumToWriteAtOnce;
					io.seek(posEntries + iostate.writeIndex*entryBytes);
					io.write(writeBuf.array(), 0, numWrittenPass1*entryBytes);
					iostate.advanceWrite(numWrittenPass1);

					int numWrittenPass2 = writtenEntries - (int)maxNumToWriteAtOnce;
					io.seek(posEntries + iostate.writeIndex*entryBytes);
					io.write(writeBuf.array(), numWrittenPass1*entryBytes, numWrittenPass2*entryBytes);
					iostate.advanceWrite(numWrittenPass2);
				}
				writeBuf.clear();
				writtenEntries = 0;

				// write zSumMax
				io.seek(posZStats + BigExp.NumBytes*confSpace.states.size());
				for (MultiStateConfSpace.State state : confSpace.states) {
					iostate.writeZSumMax[state.index].writeTo(io);
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
	 * Finishes this step and prepares for the next step.
	 * All writes are flushed to the underlying storage device by the time this method returns.
	 */
	public void finishStep() {

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
				iostate.readZSumMax[state.index].set(iostate.writeZSumMax[state.index]);
				iostate.readZSumMax[state.index].writeTo(io);
			}
			for (MultiStateConfSpace.State state : confSpace.states) {
				iostate.writeZSumMax[state.index].set(Double.NaN);
				iostate.writeZSumMax[state.index].writeTo(io);
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
