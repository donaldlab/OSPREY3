package edu.duke.cs.osprey.energy.compiled;


import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import org.joml.Vector3d;
import org.tomlj.TomlPosition;
import org.tomlj.TomlTable;


public interface EnergyCalculator {

	/** get the id of this energy calculator, that matches the forcefield ids in the conf space */
	String id();

	/** get the index of this forcefield in the conf space */
	int ffi();

	/** read runtime settings from the TOML table */
	default void readSettings(TomlTable toml, TomlPosition pos) {}

	/** calculate position-pair energy */
	double calcEnergy(double r, double r2, double[] params);

	/** get the internal energy of the static atoms */
	default double getEnergyStatic(ConfSpace.AssignedCoords coords) {
		int ffi = ffi();
		return coords.getStaticEnergy(ffi);
	}

	/** calculate the single energy of position i */
	default double calcEnergySingle(ConfSpace.AssignedCoords coords, int posi) {

		int ffi = ffi();
		ConfSpace.IndicesSingle indices = coords.getIndices(ffi, posi);

		// start with the internal energy
		double energy = indices.energy;

		// TODO: hopefully escape analysis will allocte this on the stack?
		Vector3d pos1 = new Vector3d();
		Vector3d pos2 = new Vector3d();

		// add the internal interactions
		for (int i=0; i<indices.sizeInternals(); i++) {
			int confAtom1i = indices.getInternalConfAtom1Index(i);
			int confAtom2i = indices.getInternalConfAtom2Index(i);
			int paramsi = indices.getInternalParamsIndex(i);
			coords.getConfCoords(posi, confAtom1i, pos1);
			coords.getConfCoords(posi, confAtom2i, pos2);
			double r2 = pos1.distanceSquared(pos2);
			double r = Math.sqrt(r2);
			energy += calcEnergy(r, r2, coords.getParams(ffi, paramsi));
		}

		return energy;
	}

	/** calculate the pair energy between position i and the static atoms */
	default double calcEnergyStatic(ConfSpace.AssignedCoords coords, int posi) {

		double energy = 0.0;

		// TODO: hopefully escape analysis will allocte this on the stack?
		Vector3d pos1 = new Vector3d();
		Vector3d pos2 = new Vector3d();

		int ffi = ffi();
		ConfSpace.IndicesSingle indices = coords.getIndices(ffi, posi);
		for (int i=0; i<indices.sizeStatics(); i++) {
			int confAtomi = indices.getStaticConfAtomIndex(i);
			int staticAtomi = indices.getStaticStaticAtomIndex(i);
			int paramsi = indices.getStaticParamsIndex(i);
			coords.getConfCoords(posi, confAtomi, pos1);
			coords.getStaticCoords(staticAtomi, pos2);
			double r2 = pos1.distanceSquared(pos2);
			double r = Math.sqrt(r2);
			energy += calcEnergy(r, r2, coords.getParams(ffi, paramsi));
		}

		return energy;
	}

	/** calculate the pair energy between position i1 and position i2 */
	default double calcEnergyPair(ConfSpace.AssignedCoords coords, int posi1, int posi2) {

		double energy = 0.0;

		// TODO: hopefully escape analysis will allocte this on the stack?
		Vector3d pos1 = new Vector3d();
		Vector3d pos2 = new Vector3d();

		int ffi = ffi();
		ConfSpace.IndicesPair indices = coords.getIndices(ffi, posi1, posi2);
		for (int i=0; i<indices.size(); i++) {
			int confAtomi1 = indices.getConfAtom1Index(i);
			int confAtomi2 = indices.getConfAtom2Index(i);
			int paramsi = indices.getParamsIndex(i);
			coords.getConfCoords(posi1, confAtomi1, pos1);
			coords.getConfCoords(posi2, confAtomi2, pos2);
			double r2 = pos1.distanceSquared(pos2);
			double r = Math.sqrt(r2);
			energy += calcEnergy(r, r2, coords.getParams(ffi, paramsi));
		}

		return energy;
	}

	default double calcEnergy(ConfSpace.AssignedCoords coords) {

		// start with the static energy
		double energy = getEnergyStatic(coords);

		// add the singles
		int numPos = coords.getConfSpace().positions.length;
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
}
