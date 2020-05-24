package org.mapdb.volume;

import org.mapdb.CC;
import org.mapdb.DBException;
import org.mapdb.DataIO;

import java.io.File;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;


/**
 * Based on ByteBufferMemoryVol, but with a fixed size
 */
public class FixedMemVolume extends ByteBufferVol {

	public final long maxSize;

	public FixedMemVolume(long initSize, long maxSize) {
		super(false, CC.PAGE_SHIFT, false);

		long minMaxSize = (1 << sliceShift);
		if (maxSize < minMaxSize) {
			throw new IllegalArgumentException("max size " + maxSize + " too small, must be at least " + minMaxSize);
		}

		this.maxSize = maxSize;
		if (initSize != 0) {
			ensureAvailable(initSize);
		}
	}

	@Override
	public final void ensureAvailable(long offset) {

		offset = DataIO.roundUp(offset, 1L << sliceShift);

		// check against the max size
		if (offset > maxSize) {
			throw new DBException.VolumeMaxSizeExceeded(maxSize, offset);
		}

		int slicePos = (int) (offset >>> sliceShift);

		//check for most common case, this is already mapped
		if (slicePos < slices.length) {
			return;
		}

		growLock.lock();
		try {
			//check second time
			if (slicePos <= slices.length)
				return;

			int oldSize = slices.length;
			ByteBuffer[] slices2 = Arrays.copyOf(slices, slicePos);

			for (int pos = oldSize; pos < slices2.length; pos++) {
				slices2[pos] = ByteBuffer.allocate(sliceSize).order(ByteOrder.BIG_ENDIAN);
			}

			slices = slices2;
		} catch (OutOfMemoryError e) {
			throw new DBException.OutOfMemory(e);
		} finally {
			growLock.unlock();
		}
	}


	@Override
	public void truncate(long size) {
		final int maxSize = 1 + (int) (size >>> sliceShift);
		if (maxSize == slices.length)
			return;
		if (maxSize > slices.length) {
			ensureAvailable(size);
			return;
		}
		growLock.lock();
		try {
			if (maxSize >= slices.length)
				return;
			ByteBuffer[] old = slices;
			slices = Arrays.copyOf(slices, maxSize);

			//unmap remaining buffers
			for (int i = maxSize; i < old.length; i++) {
				old[i] = null;
			}

		} finally {
			growLock.unlock();
		}
	}

	@Override
	public void close() {

		if (!closed.compareAndSet(false,true))
			return;

		growLock.lock();
		try {
			Arrays.fill(slices, null);
			slices = null;
		} finally {
			growLock.unlock();
		}
	}

	@Override
	public void sync() {
		// nothing to do
	}

	@Override
	public long length() {
		return ((long)slices.length)*sliceSize;
	}

	@Override
	public boolean isReadOnly() {
		return false;
	}

	@Override
	public File getFile() {
		return null;
	}

	@Override
	public boolean getFileLocked() {
		return false;
	}
}
