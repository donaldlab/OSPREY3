package edu.duke.cs.osprey.energy.compiled;


import edu.duke.cs.osprey.tools.HashCalculator;
import org.joml.Vector3d;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

public class ForcefieldDebugger {

	public static final ForcefieldDebugger instance = new ForcefieldDebugger();

	private static class PairKey implements Comparable<PairKey> {

		final String a;
		final String b;

		PairKey(String a, String b) {
			this.a = a;
			this.b = b;
		}

		String min() {
			return a.compareTo(b) <= 0 ? a : b;
		}

		String max() {
			return a.compareTo(b) <= 0 ? b : a;
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineHashesCommutative(
				a.hashCode(),
				b.hashCode()
			);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof PairKey && equals((PairKey)other);
		}

		public boolean equals(PairKey other) {
			return (this.a.equals(other.a) && this.b.equals(other.b))
				|| (this.a.equals(other.b) && this.b.equals(other.a));
		}

		@Override
		public String toString() {
			return String.format("%s-%s", min(), max());
		}

		@Override
		public int compareTo(PairKey other) {
			return toString().compareTo(other.toString());
		}
	}

	private static class Internal {

		final PairKey key;
		final Map<String,Double> energies = new HashMap<>();

		Internal(PairKey key) {
			this.key = key;
		}

		void add(String type, double energy) {
			if (energies.containsKey(type)) {
				throw new IllegalArgumentException(String.format(
					"already have internal energy for %s, %s", key, type
				));
			}
			energies.put(type, energy);
		}
	}

	private Map<PairKey,Internal> internals = new HashMap<>();

	public void addInternal(String idA, String idB, String type, double energy) {
		internals
			.computeIfAbsent(new PairKey(idA, idB), (key) -> new Internal(key))
			.add(type, energy);
	}

	private static class Interaction {

		final AtomPair atomPair;
		final Map<String,Double> energies = new HashMap<>();

		Interaction(AtomPair atomPair) {
			this.atomPair = atomPair;
		}

		void add(String type, double energy) {
			if (energies.containsKey(type)) {
				throw new IllegalArgumentException(String.format(
					"already have interaction energy for %s, %s = %f, can't add %f",
					atomPair.key, type, energies.get(type), energy
				));
			}
			energies.put(type, energy);
		}
	}

	private Map<PairKey,Interaction> interactions = new HashMap<>();

	public static class AtomPair {

		public final PairKey key;
		public final Vector3d a;
		public final Vector3d b;
		public final double r;

		public AtomPair(
			String idA, String idB,
			double x1, double y1, double z1,
			double x2, double y2, double z2,
			double r
		) {
			this(idA, idB, new Vector3d(x1, y1, z1), new Vector3d(x2, y2, z2), r);
		}

		public AtomPair(String idA, String idB, Vector3d a, Vector3d b, double r) {
			this.key = new PairKey(idA, idB);
			this.a = a;
			this.b = b;
			this.r = r;
		}
	}

	public void addInteraction(AtomPair atomPair, String type, double energy) {
		interactions
			.computeIfAbsent(atomPair.key, (key) -> new Interaction(atomPair))
			.add(type, energy);
	}

	public void clear() {
		internals.clear();
		interactions.clear();
	}

	public void dump(File file) {

		try (FileWriter out = new FileWriter(file)) {

			Consumer<String> write = val -> {
				try {
					out.write(val);
				} catch (IOException ex) {
					throw new RuntimeException(ex);
				}
			};

			write.accept(String.format("Internals: %d\n", internals.size()));
			forEachOrdered(internals, (key, internal) -> {
				if (internal.energies.values().stream().anyMatch(v -> v != 0.0)) {
					write.accept(String.format("%8s - %-8s", key.min(), key.max()));
					forEachOrdered(internal.energies, (type, energy) -> {
						write.accept(String.format(" %8s %7.4f", type, energy));
					});
					write.accept("\n");
				}
			});

			write.accept(String.format("Interactions: %d\n", interactions.size()));
			forEachOrdered(interactions, (key, inter) -> {
				if (inter.energies.values().stream().anyMatch(v -> v != 0.0)) {
					write.accept(String.format("%8s - %-8s", key.min(), key.max()));
					write.accept(String.format(" r=%10.6f", inter.atomPair.r));
					forEachOrdered(inter.energies, (type, energy) -> {
						write.accept(String.format(" %8s %7.4f", type, energy));
					});
					write.accept("\n");
				}
			});

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}

	private <K extends Comparable<K>,V> void forEachOrdered(Map<K,V> map, BiConsumer<K,V> block) {
		map.entrySet().stream()
			.sorted(Map.Entry.comparingByKey())
			.forEach(entry -> block.accept(entry.getKey(), entry.getValue()));
	}
}
