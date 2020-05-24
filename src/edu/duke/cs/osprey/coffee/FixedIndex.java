package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.tools.MapDBTools;
import org.mapdb.*;
import org.mapdb.serializer.GroupSerializer;
import org.mapdb.volume.FixedFileVolume;
import org.mapdb.volume.FixedMemVolume;
import org.mapdb.volume.Volume;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;


/**
 * An index that allows quick queries for entries with high and low keys,
 * but discards entries with low keys when it runs out of space
 *
 * This implementation is NOT thread-safe!
 */
public class FixedIndex<K extends Comparable<K>,V> {

	public final File file;
	public final long maxBytes;
	public final Consumer<V> evictionListener;

	private final DB db;
	private final BTreeMap<K,List<V>> map;
	private final StoreDirect store;

	public FixedIndex(File file, long maxBytes, GroupSerializer<K> keySerializer, GroupSerializer<V> valueSerializer, Consumer<V> evictionListener) {

		// MapDB requires at least 2 MiB
		if (maxBytes < 2*1024*1024) {
			throw new IllegalArgumentException("NodeIndex must have at least 2 MiB of space");
		}

		this.file = file;
		this.maxBytes = maxBytes;
		this.evictionListener = evictionListener;

		// find out how many pages fit in the space
		long numPages = maxBytes/CC.PAGE_SIZE;
		long volSize = numPages*CC.PAGE_SIZE;

		Volume vol;
		if (file != null) {
			vol = new FixedFileVolume(file, volSize, volSize);
		} else {
			vol = new FixedMemVolume(volSize, volSize);
		}

		// TODO: allow re-opening existing database?

		db = DBMaker
			.volumeDB(vol, false)
			.make();
		map = db.treeMap("index")
			.keySerializer(keySerializer)
			.valueSerializer(new MapDBTools.ValuesSerializer<>(valueSerializer))
			.counterEnable()
			.create();
		store = (StoreDirect)db.getStore();
	}

	public long size() {
		return map.sizeLong();
	}

	public K lowestKey() {
		return map.firstKey();
	}

	public K highestKey() {
		return map.lastKey();
	}

	public void add(K key, V value) {

		List<V> values = null;
		boolean wasRemoved = false;

		try {

			// try to remove the values for this key
			// this can actually cause allocations and OOM exceptions
			values = map.remove(key);
			wasRemoved = true;
			values = appendValue(values, value);

			// try to add the score if there's space
			map.put(key, values);

		} catch (DBException.VolumeMaxSizeExceeded ex) {

			// out of space!
			freeUpSpace();

			try {

				// try again
				if (!wasRemoved) {
					values = map.remove(key);
					values = appendValue(values, value);
				}
				map.put(key, values);

			} catch (DBException.VolumeMaxSizeExceeded ex2) {
				throw new Error("Can't add entry to index even after freeing up space. This is a bug.");
			}
		}
	}

	private void freeUpSpace() {
		try {

			// compact the tree so we can free up some space
			// since when we try to remove, sometimes MapDB still asks for more space
			store.compact();

			// remove 10% of the worst entries to make space
			long toRemove = (long)(map.sizeLong()*0.1);
			for (int i=0; i<toRemove; i++) {
				List<V> values = map.remove(map.firstKey());
				if (evictionListener != null) {
					for (V value : values) {
						evictionListener.accept(value);
					}
				}
			}

		} catch (DBException.VolumeMaxSizeExceeded ex) {
			// this shouldn't happen, right?
			throw new Error("Paradoxically, we ran out of space while trying to free up space in the index", ex);
		}
	}

	private List<V> appendValue(List<V> values, V value) {

		if (values == null) {
			return Collections.singletonList(value);
		}

		// in general, values ArrayList intances will be exactly sized to their contents,
		// so calling add() will usually allocate more memory than we really need

		var newValues = new ArrayList<V>(values.size() + 1);
		newValues.addAll(values);
		newValues.add(value);
		return newValues;
	}

	public V remove(K key) {

		boolean wasRemoved = false;
		V value = null;
		List<V> values = null;

		try {

			// try to remove it
			values = map.remove(key);
			wasRemoved = true;
			value = pollValue(values);
			if (values != null && !values.isEmpty()) {
				map.put(key, values);
			}
			return value;

		} catch (DBException.VolumeMaxSizeExceeded ex) {

			// out of space!
			freeUpSpace();

			try {

				if (!wasRemoved) {

					// try to remove it again
					values = map.remove(key);
					value = pollValue(values);
					if (values != null && !values.isEmpty()) {
						map.put(key, values);
					}

				} else {

					assert (values != null);
					assert (value != null);

					// just need to try to put the rest of the values back
					map.put(key, values);
				}

				return value;

			} catch (DBException.VolumeMaxSizeExceeded exAgain) {
				throw new Error("Can't remove entry from index even after freeing up space. This is a bug.");
			}
		}
	}

	private V pollValue(List<V> values) {
		if (values == null) {
			return null;
		} else {
			V value = values.get(0);
			values.remove(0);
			return value;
		}
	}
}
