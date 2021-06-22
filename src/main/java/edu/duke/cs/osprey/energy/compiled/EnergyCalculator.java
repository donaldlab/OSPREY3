package edu.duke.cs.osprey.energy.compiled;


import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import org.joml.Vector3d;

import java.io.DataInput;
import java.io.IOException;
import java.util.Arrays;


/**
 * Calculates the energy of atomic coordinates according to a forcefield definition.
 */
public interface EnergyCalculator {

	enum Type {

		// NOTE: these ids are defined by the conf space compiler in the GUI code
		Amber("amber"),
		EEF1("eef1");

		public final String id;

		Type(String id) {
			this.id = id;
		}

		public static Type get(String id) {
			return Arrays.stream(values())
				.filter(type -> type.id.equals(id))
				.findAny()
				.orElseThrow(() -> new UnsupportedOperationException("forcefield type '" + id + "' is not supported"));
		}
	}


	/** get the id of this energy calculator, that matches the forcefield ids in the conf space */
	String id();

	/** get the type of this energy calculator */
	Type type();

	/** get the index of this forcefield in the conf space */
	int ffi();

	/** read runtime settings from the stream */
	default void readSettings(DataInput in) throws IOException {}

	/** calculate position-pair energy */
	double calcEnergy(double r, double r2, double[] params);

	/** get the internal energy of the static atoms */
	default double calcEnergyStatic(AssignedCoords coords) {

		int ffi = ffi();

		// start with the static energy
		double energy = coords.getStaticEnergy(ffi);

		// TODO: hopefully escape analysis will allocate this on the stack?
		Vector3d pos1 = new Vector3d();
		Vector3d pos2 = new Vector3d();

		// add the internal interactions
		ConfSpace.IndicesStatic indices = coords.getIndices(ffi);
		for (int i=0; i<indices.size(); i++) {
			int atomi1 = indices.getStaticAtom1Index(i);
			int atomi2 = indices.getStaticAtom2Index(i);
			int paramsi = indices.getParamsIndex(i);
			coords.coords.get(coords.getStaticIndex(atomi1), pos1);
			coords.coords.get(coords.getStaticIndex(atomi2), pos2);
			double r2 = pos1.distanceSquared(pos2);
			double r = Math.sqrt(r2);

			/* DEBUG
			var res1 = coords.confSpace.resInfos[coords.confSpace.staticResInfoIndices[atomi1]];
			var res2 = coords.confSpace.resInfos[coords.confSpace.staticResInfoIndices[atomi1]];
			var atomPair = new ForcefieldDebugger.AtomPair(
				res1.chainId + res1.resId + ":" + coords.confSpace.staticNames[atomi1],
				res2.chainId + res2.resId + ":" + coords.confSpace.staticNames[atomi2]
			);
			ForcefieldDebugger.instance.addInteractionCoords(
				atomPair,
				"compiled " + type().id,
				new ForcefieldDebugger.Coords(pos1, pos2, r)
			);
			ForcefieldDebugger.instance.addInteractionEnergy(
				atomPair,
				"compiled " + type().id,
				calcEnergy(r, r2, coords.getParams(ffi, paramsi))
			);
			*/

			energy += calcEnergy(r, r2, coords.getParams(ffi, paramsi));
		}

		return energy;
	}

	/** calculate the single energy of position i */
	default double calcEnergySingle(AssignedCoords coords, int posi) {

		int ffi = ffi();

		// start with the internal energy
		double energy = coords.getInternalEnergy(ffi, posi);

		// TODO: hopefully escape analysis will allocate this on the stack?
		Vector3d pos1 = new Vector3d();
		Vector3d pos2 = new Vector3d();

		// add the internal interactions
		ConfSpace.IndicesSingle indices = coords.getIndices(ffi, posi);
		for (int i=0; i<indices.sizeInternals(); i++) {
			int confAtom1i = indices.getInternalConfAtom1Index(i);
			int confAtom2i = indices.getInternalConfAtom2Index(i);
			int paramsi = indices.getInternalParamsIndex(i);
			coords.coords.get(coords.getConfIndex(posi, confAtom1i), pos1);
			coords.coords.get(coords.getConfIndex(posi, confAtom2i), pos2);
			double r2 = pos1.distanceSquared(pos2);
			double r = Math.sqrt(r2);

			/* DEBUG
			var conf = coords.confSpace.positions[posi].confs[coords.assignments[posi]];
			var res1 = coords.confSpace.resInfos[conf.atomResInfoIndices[confAtom1i]];
			var res2 = coords.confSpace.resInfos[conf.atomResInfoIndices[confAtom2i]];
			var atomPair = new ForcefieldDebugger.AtomPair(
				res1.chainId + res1.resId + ":" + conf.atomNames[confAtom1i],
				res2.chainId + res2.resId + ":" + conf.atomNames[confAtom2i]
			);
			ForcefieldDebugger.instance.addInteractionCoords(
				atomPair,
				"compiled " + type().id,
				new ForcefieldDebugger.Coords(pos1, pos2, r)
			);
			ForcefieldDebugger.instance.addInteractionEnergy(
				atomPair,
				"compiled " + type().id,
				calcEnergy(r, r2, coords.getParams(ffi, paramsi))
			);
			*/

			energy += calcEnergy(r, r2, coords.getParams(ffi, paramsi));
		}

		return energy;
	}

	/** calculate the pair energy between position i and the static atoms */
	default double calcEnergyStatic(AssignedCoords coords, int posi) {

		double energy = 0.0;

		// TODO: hopefully escape analysis will allocate this on the stack?
		Vector3d pos1 = new Vector3d();
		Vector3d pos2 = new Vector3d();

		int ffi = ffi();
		ConfSpace.IndicesSingle indices = coords.getIndices(ffi, posi);
		for (int i=0; i<indices.sizeStatics(); i++) {
			int confAtomi = indices.getStaticConfAtomIndex(i);
			int staticAtomi = indices.getStaticStaticAtomIndex(i);
			int paramsi = indices.getStaticParamsIndex(i);
			coords.coords.get(coords.getConfIndex(posi, confAtomi), pos1);
			coords.coords.get(coords.getStaticIndex(staticAtomi), pos2);
			double r2 = pos1.distanceSquared(pos2);
			double r = Math.sqrt(r2);

			/* DEBUG
			var conf = coords.confSpace.positions[posi].confs[coords.assignments[posi]];
			var res1 = coords.confSpace.resInfos[conf.atomResInfoIndices[confAtomi]];
			var res2 = coords.confSpace.resInfos[coords.confSpace.staticResInfoIndices[staticAtomi]];
			var atomPair = new ForcefieldDebugger.AtomPair(
				res1.chainId + res1.resId + ":" + conf.atomNames[confAtomi],
				res2.chainId + res2.resId + ":" + coords.confSpace.staticNames[staticAtomi]
			);
			ForcefieldDebugger.instance.addInteractionCoords(
				atomPair,
				"compiled " + type().id,
				new ForcefieldDebugger.Coords(pos1, pos2, r)
			);
			ForcefieldDebugger.instance.addInteractionEnergy(
				atomPair,
				"compiled " + type().id,
				calcEnergy(r, r2, coords.getParams(ffi, paramsi))
			);
			*/

			energy += calcEnergy(r, r2, coords.getParams(ffi, paramsi));
		}

		return energy;
	}

	/** calculate the pair energy between position i1 and position i2 */
	default double calcEnergyPair(AssignedCoords coords, int posi1, int posi2) {

		double energy = 0.0;

		// TODO: hopefully escape analysis will allocate this on the stack?
		Vector3d pos1 = new Vector3d();
		Vector3d pos2 = new Vector3d();

		int ffi = ffi();
		ConfSpace.IndicesPair indices = coords.getIndices(ffi, posi1, posi2);
		for (int i=0; i<indices.size(); i++) {
			int confAtomi1 = indices.getConfAtom1Index(i);
			int confAtomi2 = indices.getConfAtom2Index(i);
			int paramsi = indices.getParamsIndex(i);
			coords.coords.get(coords.getConfIndex(posi1, confAtomi1), pos1);
			coords.coords.get(coords.getConfIndex(posi2, confAtomi2), pos2);
			double r2 = pos1.distanceSquared(pos2);
			double r = Math.sqrt(r2);

			/* DEBUG
			var conf1 = coords.confSpace.positions[posi1].confs[coords.assignments[posi1]];
			var conf2 = coords.confSpace.positions[posi2].confs[coords.assignments[posi2]];
			var res1 = coords.confSpace.resInfos[conf1.atomResInfoIndices[confAtomi1]];
			var res2 = coords.confSpace.resInfos[conf2.atomResInfoIndices[confAtomi2]];
			var atomPair = new ForcefieldDebugger.AtomPair(
				res1.chainId + res1.resId + ":" + conf1.atomNames[confAtomi1],
				res2.chainId + res2.resId + ":" + conf2.atomNames[confAtomi2]
			);
			ForcefieldDebugger.instance.addInteractionCoords(
				atomPair,
				"compiled " + type().id,
				new ForcefieldDebugger.Coords(pos1, pos2, r)
			);
			ForcefieldDebugger.instance.addInteractionEnergy(
				atomPair,
				"compiled " + type().id,
				calcEnergy(r, r2, coords.getParams(ffi, paramsi))
			);
			*/

			energy += calcEnergy(r, r2, coords.getParams(ffi, paramsi));
		}

		return energy;
	}

	/** calculate the total energy for the conformation */
	default double calcEnergy(AssignedCoords coords) {

		// start with the static energy
		double energy = calcEnergyStatic(coords);

		// add the singles
		int numPos = coords.confSpace.positions.length;
		for (int posi=0; posi<numPos; posi++) {
			energy += calcEnergySingle(coords, posi);
			energy += calcEnergyStatic(coords, posi);
		}

		// add the pairs
		for (int posi1=0; posi1<numPos; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				energy += calcEnergyPair(coords, posi1, posi2);
			}
		}

		return energy;
	}

	/** calculate the energy of just the given position interaction */
	default double calcEnergy(AssignedCoords coords, PosInter inter) {

		double energy;
		if (inter.posi1 == inter.posi2) {
			if (inter.posi1 == PosInter.StaticPos) {

				// static energy
				energy = calcEnergyStatic(coords);

			} else {

				// pos single energy
				energy = calcEnergySingle(coords, inter.posi1);
			}
		} else if (inter.posi1 == PosInter.StaticPos) {

			// pos-static energy
			energy = calcEnergyStatic(coords, inter.posi2);

		} else if (inter.posi2 == PosInter.StaticPos) {

			// pos-static energy
			energy = calcEnergyStatic(coords, inter.posi1);

		} else {

			// pos-pos pair energy
			energy = calcEnergyPair(coords, inter.posi1, inter.posi2);
		}

		// apply weight (but not the offset)
		return inter.weight*energy;
	}
}
