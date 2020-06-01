package edu.duke.cs.osprey.coffee;


import org.mapdb.CC;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.StoreDirect;
import org.mapdb.volume.FixedFileVolume;
import org.mapdb.volume.FixedMemVolume;
import org.mapdb.volume.Volume;

import java.io.File;


/**
 * A MapDB database that uses a fixed-size volume
 */
public class FixedDB {

	public final File file;
	public final long maxBytes;

	public final DB mapdb;
	public final StoreDirect store;

	public FixedDB(File file, long maxBytes) {

		// MapDB requires at least 2 MiB
		if (maxBytes < 2*1024*1024) {
			throw new IllegalArgumentException("FixedDB must have at least 2 MiB of space");
		}

		this.file = file;
		this.maxBytes = maxBytes;

		// find out how many pages fit in the space
		long numPages = maxBytes/ CC.PAGE_SIZE;
		long volSize = numPages*CC.PAGE_SIZE;

		Volume vol;
		if (file != null) {
			vol = new FixedFileVolume(file, volSize, volSize);
		} else {
			vol = new FixedMemVolume(volSize, volSize);
		}

		mapdb = DBMaker
			.volumeDB(vol, false)
			.make();
		store = (StoreDirect) mapdb.getStore();
	}
}
