package edu.duke.cs.osprey.confspace.compiled;

import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.compiled.motions.DihedralAngle;
import edu.duke.cs.osprey.energy.compiled.AmberEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.EEF1EnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.EnergyCalculator;
import edu.duke.cs.osprey.tools.LZMA2;

import java.io.ByteArrayInputStream;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.NoSuchElementException;


/**
 * A conformation space that reads the output produced by the ConfSpaceCompiler in the GUI.
 */
public class ConfSpace implements ConfSpaceIteration {

	public static int NotAssigned = -1;


	/** a conformation at a design position */
	public class Conf {

		public final int index;
		public final String id;
		public final String type; // TODO: use integers for types? with a lookup table to the name?
		public final int fragIndex;

		public final int numAtoms;
		public final CoordsList coords;
		public final String[] atomNames;
		public final ContinuousMotion.Description[] motions;
		/** indexed by ffi */
		public final double[] energies;

		public Conf(int index, String id, String type, int fragIndex, int numAtoms, CoordsList coords, String[] atomNames, ContinuousMotion.Description[] motions, double[] energies) {
			this.index = index;
			this.id = id;
			this.type = type;
			this.fragIndex = fragIndex;
			this.numAtoms = numAtoms;
			this.coords = coords;
			this.atomNames = atomNames;
			this.motions = motions;
			this.energies = energies;
		}
	}

	/** a design position */
	public class Pos {

		public final int index;
		public final String name;
		public final int numFrags;

		public final Conf[] confs;
		public final int maxNumAtoms;

		public Pos(int index, String name, int numFrags, Conf[] confs) {

			this.index = index;
			this.name = name;
			this.numFrags = numFrags;
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

		/** indexed by param, [0=confAtom1i, 1=confAtom2i, 2=parami] */
		private final int[][] internals;

		/** indexed by param, [0=confAtomi, 1=staticAtomi, 2=parami] */
		private final int[][] statics;

		IndicesSingle(int[][] internals, int[][] statics) {
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

	/** indexed by ff, pos, frag */
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

	/** indexed by ff, pos1:pos2, frag1, frag2 */
	private final IndicesPair[][][][] indicesPairs;


	/**
	 * Stores the actual forcefield parameters.
	 * Indexed by ff, parami (where parami comes from the indices arrays), varies (internal detail to forcefield)
	 */
	private final double[][][] ffparams;


	public static ConfSpace fromBytes(byte[] bytes) {

		// is the compiled conformation space compressed?
		// look for XZ magic bytes to see if this conf space is compressed or not
		// see XZ file spec, 2.1.1.1. Header Magic Bytes:
		// https://tukaani.org/xz/xz-file-format.txt
		byte[] xzMagic = new byte[] { (byte)0xfd, '7', 'z', 'X', 'Z', 0x00 };
		if (bytes.length > 6
			&& bytes[0] == xzMagic[0]
			&& bytes[1] == xzMagic[1]
			&& bytes[2] == xzMagic[2]
			&& bytes[3] == xzMagic[3]
			&& bytes[4] == xzMagic[4]
			&& bytes[5] == xzMagic[5]
		) {

			// yup, decompres it
			bytes = LZMA2.decompressBytes(bytes);
		}

		try (DataInputStream in = new DataInputStream(new ByteArrayInputStream(bytes))) {

			// look for the compiled conf space magic bytes "_ospccs_"
			boolean startValid =
				in.read() == '_'
				&& in.read() == 'o'
				&& in.read() == 's'
				&& in.read() == 'p'
				&& in.read() == 'c'
				&& in.read() == 'c'
				&& in.read() == 's'
				&& in.read() == '_';
			if (!startValid) {
				throw new IllegalArgumentException("unrecognized compiled conformation space format");
			}

			// read the version and call the appropriate constructor
			int version = in.readInt();
			ConfSpace confSpace;
			switch (version) {
				case 1: confSpace = new ConfSpace(in); break;
				// if we need more versions in the future:
				// case 2: confSpace = new ConfSpace(in, 0); break;
				// case 3: confSpace = new ConfSpace(in, 0, 0); break;
				// etc ...
				default: throw new IllegalArgumentException("unrecognized compiled conformation space version: " + version);
			}

			// we should see more magic bytes at the end
			// if not, something got corrupted, or there's a bug
			boolean endValid =
				in.read() == '_'
					&& in.read() == 'o'
					&& in.read() == 's'
					&& in.read() == 'p'
					&& in.read() == 'c'
					&& in.read() == 'c'
					&& in.read() == 's'
					&& in.read() == '_';
			if (!endValid) {
				throw new IllegalArgumentException("compiled conformation space has been corrupted");
			}

			return confSpace;

		} catch (IOException ex) {
			throw new RuntimeException("can't read compiled conformation space", ex);
		}
	}

	/**
	 * version 1 constructor
	 */
	private ConfSpace(DataInput in)
	throws IOException {

		// read the name
		name = in.readUTF();

		// read and build the forcefield implementations
		int numForcefields = in.readInt();
		forcefieldIds = new String[numForcefields];
		ecalcs = new EnergyCalculator[numForcefields];
		for (int ffi=0; ffi<numForcefields; ffi++) {

			String id = in.readUTF();
			forcefieldIds[ffi] = id;

			// build the forcefield implementation
			// TODO: make a better way to register and instantiate ecalc implementations?
			switch (in.readUTF()) {
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
			ecalcs[ffi].readSettings(in);
		}

		// read the static coords
		numStaticAtoms = in.readInt();
		staticCoords = new CoordsList(numStaticAtoms);
		staticNames = new String[numStaticAtoms];
		for (int i=0; i<numStaticAtoms; i++) {

			// read the coords
			staticCoords.set(
				i,
				in.readDouble(),
				in.readDouble(),
				in.readDouble()
			);

			// read the name
			staticNames[i] = in.readUTF();
		}

		// read the energies for static atoms
		staticEnergies = new double[forcefieldIds.length];
		for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
			staticEnergies[ffi] = in.readDouble();
		}

		// read the design positions, if any
		int numPositions = in.readInt();
		positions = new Pos[numPositions];
		indicesSingles = new IndicesSingle[forcefieldIds.length][numPositions][];
		for (int posi=0; posi<numPositions; posi++) {

			// read the position properties
			String posName = in.readUTF();

			// read the fragments
			int numFrags = in.readInt();

			// read the conformations
			int numConfs = in.readInt();
			for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
				indicesSingles[ffi][posi] = new IndicesSingle[numConfs];
			}
			Conf[] confs = new Conf[numConfs];
			for (int confi=0; confi<numConfs; confi++) {

				// read the conf properties
				String id = in.readUTF();
				String type = in.readUTF();
				int fragIndex = in.readInt();

				// read the conformation atoms
				int numAtoms = in.readInt();
				CoordsList atomCoords = new CoordsList(numAtoms);
				String[] atomNames = new String[numAtoms];
				for(int atomi=0; atomi<numAtoms; atomi++) {
					atomCoords.set(
						atomi,
						in.readDouble(),
						in.readDouble(),
						in.readDouble()
					);
					atomNames[atomi] = in.readUTF();
				}

				// read the motions
				int numMotions = in.readInt();
				ContinuousMotion.Description[] motions = new ContinuousMotion.Description[numMotions];
				for (int motioni=0; motioni<numMotions; motioni++) {

					String motionType = in.readUTF();
					// TODO: make a better way to register and instantiate motions?
					switch (motionType) {

						// NOTE: these string values are defined by the conf space compiler in the GUI code

						case "dihedral":
							double minDegrees = in.readDouble();
							double maxDegrees = in.readDouble();
							int a = in.readInt();
							int b = in.readInt();
							int c = in.readInt();
							int d = in.readInt();
							int numRotated = in.readInt();
							int[] rotated = new int[numRotated];
							for (int roti=0; roti<numRotated; roti++) {
								rotated[roti] = in.readInt();
							}
							motions[motioni] = new DihedralAngle.Description(
								minDegrees, maxDegrees,
								a, b, c, d,
								rotated
							);
							break;

						default:
							throw new UnsupportedOperationException("continuous motion type '" + motionType + "' is not supported");
					}
				}

				// read the internal energies
				double[] energies = new double[forcefieldIds.length];
				for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
					energies[ffi] = in.readDouble();
				}

				// finally, make the conf
				confs[confi] = new Conf(
					confi,
					id,
					type,
					fragIndex,
					numAtoms,
					atomCoords,
					atomNames,
					motions,
					energies
				);
			}

			// finally, make the position
			positions[posi] = new Pos(
				posi,
				posName,
				numFrags,
				confs
			);
		}

		// read pos and pos-static forcefield params
		for (int posi=0; posi<numPositions; posi++) {
			Pos pos = positions[posi];

			for (int fragi=0; fragi<pos.numFrags; fragi++) {

				for (int ffi=0; ffi<forcefieldIds.length; ffi++) {

					// read the pos internal forcefield params
					int numSingles = in.readInt();
					int[][] singles = new int[numSingles][3];
					for (int i=0; i<numSingles; i++) {
						singles[i][0] = in.readInt(); // atomi1
						singles[i][1] = in.readInt(); // atomi2
						singles[i][2] = in.readInt(); // parami
					}

					// read the pos-static forcefield params
					int numStatics = in.readInt();
					int[][] statics = new int[numStatics][3];
					for (int i=0; i<numStatics; i++) {
						statics[i][0] = in.readInt(); // atomi
						statics[i][1] = in.readInt(); // statici
						statics[i][2] = in.readInt(); // parami
					}

					indicesSingles[ffi][posi][fragi] = new IndicesSingle(singles, statics);
				}
			}
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

				// for each fragment pair ...
				for (int fragi1=0; fragi1<pos1.numFrags; fragi1++) {
					for (int fragi2=0; fragi2<pos2.numFrags; fragi2++) {

						// for each forcefield ...
						for (int ffi=0; ffi<forcefieldIds.length; ffi++) {

							// read the atom pairs
							int numPairs = in.readInt();
							int[][] indicesPair = new int[numPairs][3];
							for (int i=0; i<numPairs; i++) {
								indicesPair[i][0] = in.readInt(); // atomi1
								indicesPair[i][1] = in.readInt(); // atomi2
								indicesPair[i][2] = in.readInt(); // parami
							}

							indicesPairs[ffi][posPairIndex][fragi1][fragi2] = new IndicesPair(indicesPair);
						}
					}
				}
			}
		}

		// finally, read the forcefield parameters themselves
		ffparams = new double[forcefieldIds.length][][];
		for (int ffi=0; ffi<forcefieldIds.length; ffi++) {

			int numParams = in.readInt();
			ffparams[ffi] = new double[numParams][];
			for (int i=0; i<numParams; i++) {

				int numNumbers = in.readInt();
				double[] params = new double[numNumbers];
				for (int p=0; p<numNumbers; p++) {
					params[p] = in.readDouble();
				}

				ffparams[ffi][i] = params;
			}
		}
	}

	/*
	 * If we need a version 2 constructor in the future, it would look like this:
	 *
	 * private ConfSpace(DataInput in, int a)
	 * throws IOException { }
	 *
	 * version 3 would look like:
	 *
	 * private ConfSpace(DataInput in, int a, int b)
	 * throws IOException { }
	 *
	 * And so on ...
	 *
	 * Java doesn't support named constructors, so we have to differentiate by the function signature instead.
	 * The extra integers are just dummy arguments to give us different constructor signatures,
	 * but those extra arguments aren't needed for anything else.
	 */

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

	public IndicesSingle indicesSingles(int ffi, int posi, int confi) {

		// convert the conf index to a frag index
		int fragi = positions[posi].confs[confi].fragIndex;

		return indicesSingles[ffi][posi][fragi];
	}

	public IndicesPair indicesPairs(int ffi, int posi1, int confi1, int posi2, int confi2) {

		// convert the conf indices to frag indices
		int fragi1 = positions[posi1].confs[confi1].fragIndex;
		int fragi2 = positions[posi2].confs[confi2].fragIndex;

		return indicesPairs[ffi][posPairIndex(posi1, posi2)][fragi1][fragi2];
	}

	public double[] ffparams(int ffi, int paramsi) {
		return ffparams[ffi][paramsi];
	}

	@Override
	public int countSingles() {
		int count = 0;
		for (int posi1=1; posi1<positions.length; posi1++) {
			count += positions[posi1].confs.length;
		}
		return count;
	}

	@Override
	public int countPairs() {
		int count = 0;
		for (int posi1=1; posi1<positions.length; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				count += positions[posi1].confs.length*positions[posi2].confs.length;
			}
		}
		return count;
	}

	@Override
	public int numPos() {
		return positions.length;
	}

	@Override
	public int numConf(int posi) {
		return positions[posi].confs.length;
	}

	@Override
	public String confResType(int posi, int confi) {
		return positions[posi].confs[confi].type;
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

	public int[] assign(RCTuple tuple) {
		int[] assignments = assign();
		for (int i=0; i<tuple.size(); i++) {
			assignments[tuple.pos.get(i)] = tuple.RCs.get(i);
		}
		return assignments;
	}

	/** Makes coordinates with the given assignments */
	public AssignedCoords makeCoords(int[] assignments) {
		return new AssignedCoords(this, assignments);
	}
}
