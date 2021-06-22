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

	public static class AtomPair implements Comparable<AtomPair> {

		final String a;
		final String b;

		public AtomPair(String a, String b) {
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
			return other instanceof AtomPair && equals((AtomPair)other);
		}

		public boolean equals(AtomPair other) {
			return (this.a.equals(other.a) && this.b.equals(other.b))
				|| (this.a.equals(other.b) && this.b.equals(other.a));
		}

		@Override
		public String toString() {
			return String.format("%s-%s", min(), max());
		}

		@Override
		public int compareTo(AtomPair other) {
			return toString().compareTo(other.toString());
		}
	}

	public static class Coords {

		public final Vector3d a;
		public final Vector3d b;
		public final double r;

		public Coords(Vector3d a, Vector3d b, double r) {
			this.a = new Vector3d(a);
			this.b = new Vector3d(b);
			this.r = r;
		}

		public Coords(
			double x1, double y1, double z1,
			double x2, double y2, double z2,
			double r
		) {
			this.a = new Vector3d(x1, y1, z1);
			this.b = new Vector3d(x2, y2, z2);
			this.r = r;
		}
	}


	private static class Internal {

		final AtomPair key;
		final Map<String,Double> energies = new HashMap<>();

		Internal(AtomPair key) {
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

	private Map<AtomPair,Internal> internals = new HashMap<>();

	public void addInternal(String idA, String idB, String type, double energy) {
		internals
			.computeIfAbsent(new AtomPair(idA, idB), (key) -> new Internal(key))
			.add(type, energy);
	}
	
	private static class Interaction {

		final AtomPair atomPair;
		final Map<String,Coords> coords = new HashMap<>();
		final Map<String,Double> energies = new HashMap<>();

		Interaction(AtomPair atomPair) {
			this.atomPair = atomPair;
		}

		void addCoords(String type, Coords coords) {
			if (this.coords.containsKey(type)) {
				throw new IllegalArgumentException(String.format(
					"already have interaction coords for %s, %s, can't add more",
					atomPair, type
				));
			}
			this.coords.put(type, coords);
		}

		void addEnergy(String type, double energy) {
			if (energies.containsKey(type)) {
				throw new IllegalArgumentException(String.format(
					"already have interaction energy for %s, %s = %f, can't add %f",
					atomPair, type, energies.get(type), energy
				));
			}
			energies.put(type, energy);
		}
	}

	private Map<AtomPair,Interaction> interactions = new HashMap<>();

	public void addInteractionCoords(AtomPair atomPair, String type, Coords coords) {
		interactions
			.computeIfAbsent(atomPair, (key) -> new Interaction(atomPair))
			.addCoords(type, coords);
	}

	public void addInteractionEnergy(AtomPair atomPair, String type, double energy) {
		interactions
			.computeIfAbsent(atomPair, (key) -> new Interaction(atomPair))
			.addEnergy(type, energy);
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
					forEachOrdered(inter.coords, (type, coords) -> {
						write.accept(String.format(" %8s (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f) r=%10.6f",
							type,
							coords.a.x, coords.a.y, coords.a.z,
							coords.b.x, coords.b.y, coords.b.z,
							coords.r
						));
					});
					forEachOrdered(inter.energies, (type, energy) -> {
						write.accept(String.format(" %8s e=%7.4f", type, energy));
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
