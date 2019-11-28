package edu.duke.cs.osprey.confspace.compiled;


import edu.duke.cs.osprey.confspace.compiled.motions.DihedralAngle;
import edu.duke.cs.osprey.energy.compiled.AmberEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.EEF1EnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.EnergyCalculator;
import edu.duke.cs.osprey.tools.TomlTools;
import org.tomlj.*;

import java.util.Arrays;
import java.util.NoSuchElementException;


/**
 * A conformation space that reads the output produced by the ConfSpaceCompiler in the GUI.
 */
public class ConfSpace {

	public static int NotAssigned = -1;


	/** a conformation at a design position */
	public class Conf {

		public final int index;
		public final String id;
		public final String type; // TODO: use integers for types? with a lookup table to the name?

		public final int numAtoms;
		public final CoordsList coords;
		public final String[] atomNames;

		public final ContinuousMotion.Description[] motions;

		public Conf(int index, String id, String type, int numAtoms, CoordsList coords, String[] atomNames, ContinuousMotion.Description[] motions) {
			this.index = index;
			this.id = id;
			this.type = type;
			this.numAtoms = numAtoms;
			this.coords = coords;
			this.atomNames = atomNames;
			this.motions = motions;
		}
	}

	/** a design position */
	public class Pos {

		public final int index;
		public final String name;

		public final Conf[] confs;
		public final int maxNumAtoms;

		public Pos(int index, String name, Conf[] confs) {

			this.index = index;
			this.name = name;
			this.confs = confs;

			maxNumAtoms = Arrays.stream(confs)
				.mapToInt(conf -> conf.numAtoms)
				.max()
				.orElse(0);
		}

		public Conf findConf(String id) {
			for (Conf conf : confs) {
				if (conf.id.equals(id)) {
					return conf;
				}
			}
			return null;
		}

		public Conf findConfOrThrow(String id) {
			Conf conf = findConf(id);
			if (conf != null) {
				return conf;
			}
			throw new NoSuchElementException("no conf with id " + id);
		}
	}

	public final String name;
	public final String[] forcefieldIds;
	public final EnergyCalculator[] ecalcs;

	public final int numStaticAtoms;
	public final CoordsList staticCoords;
	public final String[] staticNames;
	public final double[] staticEnergies;

	public final Pos[] positions;

	public class IndicesSingle {

		public final double energy;

		/** indexed by param, [0=confAtom1i, 1=confAtom2i, 2=parami] */
		private final int[][] internals;

		/** indexed by param, [0=confAtomi, 1=staticAtomi, 2=parami] */
		private final int[][] statics;

		IndicesSingle(double energy, int[][] internals, int[][] statics) {
			this.energy = energy;
			this.internals = internals;
			this.statics = statics;
		}

		public int sizeInternals() {
			return internals.length;
		}
		public int getInternalConfAtom1Index(int i) {
			return internals[i][0];
		}
		public int getInternalConfAtom2Index(int i) {
			return internals[i][1];
		}
		public int getInternalParamsIndex(int i) {
			return internals[i][2];
		}

		public int sizeStatics() {
			return statics.length;
		}
		public int getStaticConfAtomIndex(int i) {
			return statics[i][0];
		}
		public int getStaticStaticAtomIndex(int i) {
			return statics[i][1];
		}
		public int getStaticParamsIndex(int i) {
			return statics[i][2];
		}
	}

	/** indexed by ff, pos, conf */
	private final IndicesSingle[][][] indicesSingles;

	public class IndicesPair {

		/** indexed by param, [0=conf1Atomi, 1=conf2Atomi, 2=parami] */
		private final int[][] indices;

		public IndicesPair(int[][] indices) {
			this.indices = indices;
		}

		public int size() {
			return indices.length;
		}
		public int getConfAtom1Index(int i) {
			return indices[i][0];
		}
		public int getConfAtom2Index(int i) {
			return indices[i][1];
		}
		public int getParamsIndex(int i) {
			return indices[i][2];
		}
	}

	/** indexed by ff, pos1:pos2, conf1, conf2 */
	private final IndicesPair[][][][] indicesPairs;


	/**
	 * Stores the actual forcefield parameters.
	 * Indexed by ff, parami (where parami comes from the indices arrays), varies (internal detail to forcefield)
	 */
	private final double[][][] ffparams;

	public ConfSpace(String toml) {

		// parse the TOML file
		TomlTable doc = TomlTools.parseOrThrow(toml);

		name = TomlTools.getStringOrThrow(doc, "name");

		// read and build the forcefield implementations
		TomlArray forcefieldsArray = TomlTools.getArrayOrThrow(doc, "forcefields");
		forcefieldIds = new String[forcefieldsArray.size()];
		ecalcs = new EnergyCalculator[forcefieldsArray.size()];
		for (int ffi=0; ffi<forcefieldsArray.size(); ffi++) {
			TomlArray ffArray = TomlTools.getArrayOrThrow(forcefieldsArray, ffi);
			TomlPosition ffPos = forcefieldsArray.inputPositionOf(ffi);

			String id = TomlTools.getStringOrThrow(ffArray, 0, ffPos);
			forcefieldIds[ffi] = id;

			// TODO: make a better way to register and instantiate ecalc implementations?
			String implementation = TomlTools.getStringOrThrow(ffArray, 1, ffPos);
			switch (implementation) {
				case "amber":
					ecalcs[ffi] = new AmberEnergyCalculator(id, ffi);
				break;
				case "eef1":
					ecalcs[ffi] = new EEF1EnergyCalculator(id, ffi);
				break;
			}
		}

		// read the forcefield settings, when needed
		for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
			String ffkey = "ffsettings." + ffi;
			TomlTable ffsettingsTable = doc.getTable(ffkey);
			TomlPosition ffsettingsPos = doc.inputPositionOf(ffkey);
			if (ffsettingsTable != null) {
				ecalcs[ffi].readSettings(ffsettingsTable, ffsettingsPos);
			}
		}

		// read the static coords
		TomlTable staticTable = TomlTools.getTableOrThrow(doc, "static");
		TomlArray staticAtomsArray = TomlTools.getArrayOrThrow(staticTable, "atoms");
		numStaticAtoms = staticAtomsArray.size();
		staticCoords = new CoordsList(numStaticAtoms);
		staticNames = new String[numStaticAtoms];
		for (int i=0; i<numStaticAtoms; i++) {
			TomlTable table = TomlTools.getTableOrThrow(staticAtomsArray, i);
			TomlPosition tablePos = staticAtomsArray.inputPositionOf(i);

			// read the coords
			TomlArray coordsArray = TomlTools.getArrayOrThrow(table, "xyz", tablePos);
			staticCoords.set(
				i,
				coordsArray.getDouble(0),
				coordsArray.getDouble(1),
				coordsArray.getDouble(2)
			);

			// read the name
			staticNames[i] = TomlTools.getStringOrThrow(table, "name", tablePos);
		}

		// read the static energies
		TomlTable staticEnergyTable = TomlTools.getTableOrThrow(staticTable, "energy");
		TomlPosition staticEnergyPos = staticTable.inputPositionOf("energy");
		staticEnergies = new double[forcefieldIds.length];
		for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
			String ffkey = "" + ffi;
			staticEnergies[ffi] = TomlTools.getDoubleOrThrow(staticEnergyTable, ffkey, staticEnergyPos);
		}

		// read the design positions, if any
		TomlTable positionsTable = doc.getTable("pos");
		int numPositions;
		if (positionsTable != null) {
			numPositions = positionsTable.size();
		} else {
			numPositions = 0;
		}

		positions = new Pos[numPositions];
		indicesSingles = new IndicesSingle[forcefieldIds.length][numPositions][];

		for (int posi=0; posi<numPositions; posi++) {
			String posKey = "" + posi;
			TomlTable posTable = TomlTools.getTableOrThrow(positionsTable, posKey);
			TomlPosition posPos = positionsTable.inputPositionOf(posKey);

			String posName = TomlTools.getStringOrThrow(posTable, "name", posPos);

			// read the conformations
			TomlTable confsTable = TomlTools.getTableOrThrow(posTable, "conf", posPos);
			int numConfs = confsTable.size();

			for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
				indicesSingles[ffi][posi] = new IndicesSingle[numConfs];
			}

			Conf[] confs = new Conf[numConfs];
			for (int confi=0; confi<numConfs; confi++) {
				String confKey = "" + confi;
				TomlTable confTable = TomlTools.getTableOrThrow(confsTable, confKey, posPos);
				TomlPosition confPos = confsTable.inputPositionOf(confKey);

				String id = TomlTools.getStringOrThrow(confTable, "id", confPos);
				String type = TomlTools.getStringOrThrow(confTable, "type", confPos);

				// read the conformation atoms
				TomlArray atomsArray = TomlTools.getArrayOrThrow(confTable, "atoms", confPos);
				int numAtoms = atomsArray.size();
				CoordsList atomCoords = new CoordsList(numAtoms);
				String[] atomNames = new String[numAtoms];
				for(int atomi=0; atomi<numAtoms; atomi++) {
					TomlTable atomTable = TomlTools.getTableOrThrow(atomsArray, atomi);
					TomlPosition atomPos = atomsArray.inputPositionOf(atomi);

					atomNames[atomi] = TomlTools.getStringOrThrow(atomTable, "name", atomPos);

					// read the coords
					TomlArray xyzArray = TomlTools.getArrayOrThrow(atomTable, "xyz", atomPos);
					atomCoords.set(
						atomi,
						TomlTools.getDoubleOrThrow(xyzArray, 0, atomPos),
						TomlTools.getDoubleOrThrow(xyzArray, 1, atomPos),
						TomlTools.getDoubleOrThrow(xyzArray, 2, atomPos)
					);
				}

				// read the DoF descriptions
				ContinuousMotion.Description[] motions;
				TomlArray motionsArray = confTable.getArray("motions");
				if (motionsArray == null) {

					motions = new ContinuousMotion.Description[0];

				} else {

					TomlPosition motionsPos = confTable.inputPositionOf("motions");
					motions = new ContinuousMotion.Description[motionsArray.size()];
					for (int i=0; i<motionsArray.size(); i++) {
						TomlTable motionTable = TomlTools.getTableOrThrow(motionsArray, i, motionsPos);
						TomlPosition motionPos = motionsArray.inputPositionOf(i);

						// TODO: make a better way to register and instantiate motions?
						String motionType = TomlTools.getStringOrThrow(motionTable, "type", motionPos);
						switch (motionType) {

							case "dihedral":
								TomlArray boundsArray = TomlTools.getArrayOrThrow(motionTable, "bounds", motionPos);
								TomlArray abcdArray = TomlTools.getArrayOrThrow(motionTable, "abcd", motionPos);
								TomlArray rotatedArray = TomlTools.getArrayOrThrow(motionTable, "rotated", motionPos);
								int[] rotated = new int[rotatedArray.size()];
								for (int j=0; j<rotatedArray.size(); j++) {
									rotated[j] = TomlTools.getIntOrThrow(rotatedArray, j, motionPos);
								}
								motions[i] = new DihedralAngle.Description(
									TomlTools.getDoubleOrThrow(boundsArray, 0, motionPos),
									TomlTools.getDoubleOrThrow(boundsArray, 1, motionPos),
									TomlTools.getIntOrThrow(abcdArray, 0, motionPos),
									TomlTools.getIntOrThrow(abcdArray, 1, motionPos),
									TomlTools.getIntOrThrow(abcdArray, 2, motionPos),
									TomlTools.getIntOrThrow(abcdArray, 3, motionPos),
									rotated
								);
							break;

							default:
								throw new UnsupportedOperationException("continuous motion type '" + motionType + "' is not supported");
						}
					}
				}

				// finally, make the conf
				confs[confi] = new Conf(
					confi,
					id,
					type,
					numAtoms,
					atomCoords,
					atomNames,
					motions
				);

				// read the internal energies
				TomlTable energyTable = TomlTools.getTableOrThrow(confTable, "energy", confPos);
				TomlPosition energyPos = confTable.inputPositionOf("energy");

				// read the forcefield params
				TomlTable paramsInternalTable = TomlTools.getTableOrThrow(confTable, "params.internal", confPos);
				TomlPosition paramsInternalPos = confTable.inputPositionOf("params.internal");
				TomlTable paramsStaticTable = TomlTools.getTableOrThrow(confTable, "params.static", confPos);
				TomlPosition paramsStaticPos = confTable.inputPositionOf("params.static");

				for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
					String ffkey = "" + ffi;

					// read the energy
					double energy = TomlTools.getDoubleOrThrow(energyTable, ffkey, energyPos);

					// read the pos internal forcefield params
					TomlArray internalArray = TomlTools.getArrayOrThrow(paramsInternalTable, ffkey, paramsInternalPos);
					TomlPosition internalPos = paramsInternalTable.inputPositionOf(ffkey);
					int[][] singles = new int[internalArray.size()][3];
					for (int i=0; i<internalArray.size(); i++) {
						TomlArray indicesArray = TomlTools.getArrayOrThrow(internalArray, i, internalPos);
						TomlPosition indicesPos = internalArray.inputPositionOf(i);
						singles[i][0] = TomlTools.getIntOrThrow(indicesArray, 0, indicesPos); // atom1i
						singles[i][1] = TomlTools.getIntOrThrow(indicesArray, 1, indicesPos); // atom2i
						singles[i][2] = TomlTools.getIntOrThrow(indicesArray, 2, indicesPos); // parami
					}

					// read the pos-static forcefield params
					TomlArray staticArray = TomlTools.getArrayOrThrow(paramsStaticTable, ffkey, paramsStaticPos);
					TomlPosition staticPos = paramsStaticTable.inputPositionOf(ffkey);
					int[][] statics = new int[staticArray.size()][3];
					for (int i=0; i<staticArray.size(); i++) {
						TomlArray indicesArray = TomlTools.getArrayOrThrow(staticArray, i, staticPos);
						TomlPosition indicesPos = staticArray.inputPositionOf(i);
						statics[i][0] = TomlTools.getIntOrThrow(indicesArray, 0, indicesPos); // atomi
						statics[i][1] = TomlTools.getIntOrThrow(indicesArray, 1, indicesPos); // statici
						statics[i][2] = TomlTools.getIntOrThrow(indicesArray, 2, indicesPos); // parami
					}

					indicesSingles[ffi][posi][confi] = new IndicesSingle(energy, singles, statics);
				}
			}

			// finally, make the position
			positions[posi] = new Pos(
				posi,
				posName,
				confs
			);
		}

		// read pos-pos forcefield params
		int numPosPairs = Math.max(0, numPositions*(numPositions + 1)/2 - 1);
		indicesPairs = new IndicesPair[forcefieldIds.length][numPosPairs][][];

		// for each position pair ...
		for (int posi1=0; posi1<numPositions; posi1++) {
			Pos pos1 = positions[posi1];
			for (int posi2=0; posi2<posi1; posi2++) {
				Pos pos2 = positions[posi2];

				int posPairIndex = posPairIndex(pos1.index, pos2.index);

				for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
					indicesPairs[ffi][posPairIndex] = new IndicesPair[pos1.confs.length][pos2.confs.length];
				}

				// for each conf pair ...
				for (Conf conf1 : pos1.confs) {
					for (Conf conf2 : pos2.confs) {

						String key = String.format("pos.%d.conf.%d.params.pos.%d.conf.%d",
							pos1.index, conf1.index, pos2.index, conf2.index
						);
						TomlTable pairTable = TomlTools.getTableOrThrow(doc, key);
						TomlPosition pairPos = doc.inputPositionOf(key);

						// for each forcefield ...
						for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
							String ffkey = "" + ffi;
							TomlArray ffArray = TomlTools.getArrayOrThrow(pairTable, ffkey, pairPos);
							TomlPosition ffPos = pairTable.inputPositionOf(ffkey);

							int[][] indicesPair = new int[ffArray.size()][3];
							for (int i=0; i<ffArray.size(); i++) {
								TomlArray indicesArray = TomlTools.getArrayOrThrow(ffArray, i, ffPos);
								TomlPosition indicesPos = ffArray.inputPositionOf(i);
								indicesPair[i][0] = TomlTools.getIntOrThrow(indicesArray, 0, indicesPos); // atomi1
								indicesPair[i][1] = TomlTools.getIntOrThrow(indicesArray, 1, indicesPos); // atomi2
								indicesPair[i][2] = TomlTools.getIntOrThrow(indicesArray, 2, indicesPos); // parami
							}

							indicesPairs[ffi][posPairIndex][conf1.index][conf2.index] = new IndicesPair(indicesPair);
						}
					}
				}
			}
		}

		// finally, read the parameters themselves
		TomlTable ffparamsTable = TomlTools.getTableOrThrow(doc, "ffparams");
		TomlPosition ffparamsPos = doc.inputPositionOf("ffparams");
		ffparams = new double[forcefieldIds.length][][];
		for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
			String ffkey = "" + ffi;
			TomlArray ffArray = TomlTools.getArrayOrThrow(ffparamsTable, ffkey, ffparamsPos);
			TomlPosition ffPos = ffparamsTable.inputPositionOf(ffkey);

			int numParams = ffArray.size();
			ffparams[ffi] = new double[numParams][];
			for (int i=0; i<numParams; i++) {
				TomlArray paramsArray = TomlTools.getArrayOrThrow(ffArray, i, ffPos);
				TomlPosition paramsPos = ffArray.inputPositionOf(i);

				double[] params = new double[paramsArray.size()];
				for (int p=0; p<paramsArray.size(); p++) {
					params[p] = TomlTools.getDoubleOrThrow(paramsArray, p, paramsPos);
				}

				ffparams[ffi][i] = params;
			}
		}
	}

	private int posPairIndex(int posi1, int posi2) {

		// posi2 should be strictly less than posi1
		if (posi2 > posi1) {
			int swap = posi1;
			posi1 = posi2;
			posi2 = swap;
		} else if (posi1 == posi2) {
			throw new Error("Can't pair position " + posi1 + " with itself");
		}

		return posi1*(posi1 - 1)/2 + posi2;
	}

	public int[] getNumConfsAtPos() {
		return Arrays.stream(positions)
			.mapToInt(pos -> pos.confs.length)
			.toArray();
	}

	public IndicesSingle indicesSingles(int ffi, int posi, int confi) {
		return indicesSingles[ffi][posi][confi];
	}

	public IndicesPair indicesPairs(int ffi, int posi1, int confi1, int posi2, int confi2) {
		return indicesPairs[ffi][posPairIndex(posi1, posi2)][confi1][confi2];
	}

	public double[] ffparams(int ffi, int paramsi) {
		return ffparams[ffi][paramsi];
	}

	/** Gets the total number of confs for all positions */
	public int countConfSingles() {
		int count = 0;
		for (int posi1=1; posi1<positions.length; posi1++) {
			count += positions[posi1].confs.length;
		}
		return count;
	}

	/** Gets the total number of conf pairs for all positions */
	public int countConfPairs() {
		int count = 0;
		for (int posi1=1; posi1<positions.length; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				count += positions[posi1].confs.length*positions[posi2].confs.length;
			}
		}
		return count;
	}

	/** Makes assignments with no assigned positions. */
	public int[] assign() {
		int[] assignments = new int[positions.length];
		Arrays.fill(assignments, NotAssigned);
		return assignments;
	}

	/** Makes assignments with one assigned position */
	public int[] assign(int posi1, int confi1) {
		int[] assignments = assign();
		assignments[posi1] = confi1;
		return assignments;
	}

	/** Makes assignments with two assigned positions */
	public int[] assign(int posi1, int confi1, int posi2, int confi2) {
		int[] assignments = assign(posi1, confi1);
		assignments[posi2] = confi2;
		return assignments;
	}

	/** Makes assignments with three assigned positions */
	public int[] assign(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3) {
		int[] assignments = assign(posi1, confi1, posi2, confi2);
		assignments[posi3] = confi3;
		return assignments;
	}

	/** Makes coordinates with the given assignments */
	public AssignedCoords makeCoords(int[] assignments) {
		return new AssignedCoords(this, assignments);
	}
}
