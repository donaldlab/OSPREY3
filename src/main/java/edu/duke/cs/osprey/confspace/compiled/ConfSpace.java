package edu.duke.cs.osprey.confspace.compiled;

import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.compiled.motions.DihedralAngle;
import edu.duke.cs.osprey.confspace.compiled.motions.TranslationRotation;
import edu.duke.cs.osprey.energy.compiled.AmberEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.EEF1EnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.EnergyCalculator;
import edu.duke.cs.osprey.tools.LZMA2;
import org.joml.Vector3d;

import java.io.ByteArrayInputStream;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
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
		public final int[] atomMolInfoIndices;
		public final int[] atomResInfoIndices;
		public final ContinuousMotion.ConfDescription[] motions;
		/** indexed by ffi */
		public final double[] energies;

		public Conf(int index, String id, String type, int fragIndex, int numAtoms, CoordsList coords, String[] atomNames, int[] atomMolInfoIndices, int[] atomResInfoIndices, ContinuousMotion.ConfDescription[] motions, double[] energies) {
			this.index = index;
			this.id = id;
			this.type = type;
			this.fragIndex = fragIndex;
			this.numAtoms = numAtoms;
			this.coords = coords;
			this.atomNames = atomNames;
			this.atomResInfoIndices = atomResInfoIndices;
			this.atomMolInfoIndices = atomMolInfoIndices;
			this.motions = motions;
			this.energies = energies;
		}

		public Integer findAtomIndex(String name) {
			for (int atomi=0; atomi<numAtoms; atomi++) {
				if (atomNames[atomi].equals(name)) {
					return atomi;
				}
			}
			return null;
		}

		public int findAtomIndexOrThrow(String name) {
			Integer atomi = findAtomIndex(name);
			if (atomi != null) {
				return atomi;
			}
			throw new NoSuchElementException("no atom with name " + name);
		}
	}

	/** a design position */
	public class Pos {

		public final int index;
		public final String name;
		public final String wildType;
		public final int numFrags;

		public final Conf[] confs;
		public final int maxNumAtoms;
		public final boolean hasMutations;

		public Pos(int index, String name, String wildType, int numFrags, Conf[] confs) {

			this.index = index;
			this.name = name;
			this.wildType = wildType;
			this.numFrags = numFrags;
			this.confs = confs;

			maxNumAtoms = Arrays.stream(confs)
				.mapToInt(conf -> conf.numAtoms)
				.max()
				.orElse(0);

			hasMutations = Arrays.stream(confs)
				.anyMatch(conf -> !conf.type.equals(wildType));
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

	public static class MolInfo {

		public final String name;
		/** can be null */
		public final String type;
		public final ContinuousMotion.MolDescription[] motions;

		public MolInfo(String name, String type, int numMotions) {
			this.name = name;
			this.type = type;
			this.motions = new ContinuousMotion.MolDescription[numMotions];
		}
	}

	public static class ResInfo {

		public final String chainId;
		public final String resId;
		public final String resType;
		public final int indexInChain;

		public ResInfo(String chainId, String resId, String resType, int indexInChain) {
			this.chainId = chainId;
			this.resId = resId;
			this.resType = resType;
			this.indexInChain = indexInChain;
		}
	}

	private final int hash;

	public final String name;
	public final String[] forcefieldIds;
	public final EnergyCalculator[] ecalcs;

	public final MolInfo[] molInfos;
	public final ResInfo[] resInfos;

	public final int numStaticAtoms;
	public final CoordsList staticCoords;
	public final String[] staticNames;
	public final int[] staticMolInfoIndices;
	public final int[] staticResInfoIndices;
	public final double[] staticEnergies;

	public final Pos[] positions;
	public final int maxNumConfAtoms;
	public final int[] confAtomOffsetsByPos;
	public final int maxNumDofs;
	public final SeqSpace seqSpace;

	public class IndicesStatic {

		/** indexed by param, [0=atomi1, 1=atomi2, 2=parami] */
		private final int[][] indices;

		public IndicesStatic(int[][] indices) {
			this.indices = indices;
		}

		public int size() {
			return indices.length;
		}
		public int getStaticAtom1Index(int i) {
			return indices[i][0];
		}
		public int getStaticAtom2Index(int i) {
			return indices[i][1];
		}
		public int getParamsIndex(int i) {
			return indices[i][2];
		}
	}

	/** indexed by ff */
	private final IndicesStatic[] indicesStatic;

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

		// compute the hash code from the raw bytes
		int hash = Arrays.hashCode(bytes);

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
				case 1: confSpace = new ConfSpace(hash, in); break;
				// if we need more versions in the future:
				// case 2: confSpace = new ConfSpace(hash, in, 0); break;
				// case 3: confSpace = new ConfSpace(hash, in, 0, 0); break;
				// etc ...
				// tragically Java doesn't have named constructors, and static constructor functions don't mix well with final fields
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
	private ConfSpace(int hash, DataInput in)
	throws IOException {

		// save the hash
		this.hash = hash;

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
			ecalcs[ffi] = switch (EnergyCalculator.Type.get(in.readUTF())) {
				case Amber -> new AmberEnergyCalculator(id, ffi);
				case EEF1 -> new EEF1EnergyCalculator(id, ffi);
			};
		}

		// read the forcefield settings, when needed
		for (int ffi=0; ffi<forcefieldIds.length; ffi++) {
			ecalcs[ffi].readSettings(in);
		}

		// read the molecule infos
		molInfos = new MolInfo[in.readInt()];
		for (int i=0; i<molInfos.length; i++) {
			molInfos[i] = new MolInfo(
				in.readUTF(),
				in.readByte() == 1 ? in.readUTF() : null,
				in.readInt()
			);

			// read the motions, if any
			for (int motioni=0; motioni<molInfos[i].motions.length; motioni++) {

				String motionType = in.readUTF();
				// TODO: make a better way to register and instantiate motions?
				switch (motionType) {

					// NOTE: these string values are defined by the conf space compiler in the GUI code

					case "dihedralAngle":
						molInfos[i].motions[motioni] = readDihedralAngle(in);
					break;

					case "translationRotation":
						molInfos[i].motions[motioni] = readTranslationRotation(in);
					break;

					default:
						throw new UnsupportedOperationException("continuous motion type '" + motionType + "' is not supported for molecules");
				}
			}
		}

		// read the res infos
		resInfos = new ResInfo[in.readInt()];
		for (int i=0; i<resInfos.length; i++) {
			resInfos[i] = new ResInfo(
				in.readUTF(),
				in.readUTF(),
				in.readUTF(),
				in.readInt()
			);
		}

		// read the static coords
		numStaticAtoms = in.readInt();
		staticCoords = new CoordsList(numStaticAtoms);
		staticNames = new String[numStaticAtoms];
		staticResInfoIndices = new int[numStaticAtoms];
		staticMolInfoIndices = new int[numStaticAtoms];
		for (int i=0; i<numStaticAtoms; i++) {

			// read the coords
			staticCoords.set(
				i,
				in.readDouble(),
				in.readDouble(),
				in.readDouble()
			);

			// read the metadata
			staticMolInfoIndices[i] = in.readInt();
			staticResInfoIndices[i] = in.readInt();
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
			String wildType = in.readUTF();

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
				int[] atomMolInfoIndices = new int[numAtoms];
				int[] atomResInfoIndices = new int[numAtoms];
				for(int atomi=0; atomi<numAtoms; atomi++) {
					atomCoords.set(
						atomi,
						in.readDouble(),
						in.readDouble(),
						in.readDouble()
					);
					atomMolInfoIndices[atomi] = in.readInt();
					atomResInfoIndices[atomi] = in.readInt();
					atomNames[atomi] = in.readUTF();
				}

				// read the motions
				int numMotions = in.readInt();
				ContinuousMotion.ConfDescription[] motions = new ContinuousMotion.ConfDescription[numMotions];
				for (int motioni=0; motioni<numMotions; motioni++) {

					String motionType = in.readUTF();
					// TODO: make a better way to register and instantiate motions?
					switch (motionType) {

						// NOTE: these string values are defined by the conf space compiler in the GUI code

						case "dihedralAngle":
							motions[motioni] = readDihedralAngle(in);
						break;

						default:
							throw new UnsupportedOperationException("continuous motion type '" + motionType + "' is not supported for conformations");
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
					atomMolInfoIndices,
					atomResInfoIndices,
					motions,
					energies
				);
			}

			// finally, make the position
			positions[posi] = new Pos(
				posi,
				posName,
				wildType,
				numFrags,
				confs
			);
		}

		// compute the conf atoms offsets
		confAtomOffsetsByPos = new int[positions.length];
		int maxNumConfAtoms = staticCoords.size;
		for (ConfSpace.Pos pos : positions) {
			confAtomOffsetsByPos[pos.index] = maxNumConfAtoms;
			maxNumConfAtoms += pos.maxNumAtoms;
		}
		this.maxNumConfAtoms = maxNumConfAtoms;

		// read the static forcefield params
		indicesStatic = new IndicesStatic[forcefieldIds.length];
		for (int ffi=0; ffi<forcefieldIds.length; ffi++) {

			int num = in.readInt();
			int[][] indices = new int[num][3];
			for (int i=0; i<num; i++) {
				indices[i][0] = in.readInt(); // atomi1
				indices[i][1] = in.readInt(); // atomi2
				indices[i][2] = in.readInt(); // parami
			}

			indicesStatic[ffi] = new IndicesStatic(indices);
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

		// count the max number of dofs in a conformation
		this.maxNumDofs =
			Arrays.stream(molInfos)
				.flatMap(info -> Arrays.stream(info.motions))
				.mapToInt(motion -> motion.maxNumDofs())
				.sum()
			+ Arrays.stream(positions)
				.mapToInt(pos -> Arrays.stream(pos.confs)
					.mapToInt(conf -> Arrays.stream(conf.motions)
						.mapToInt(motion -> motion.maxNumDofs())
						.sum()
					)
					.max()
					.orElse(0)
				)
				.sum();

		seqSpace = new SeqSpace(this);
	}

	private static DihedralAngle.Description readDihedralAngle(DataInput in)
	throws IOException {
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
		return new DihedralAngle.Description(
			minDegrees, maxDegrees,
			a, b, c, d,
			rotated
		);
	}

	private static TranslationRotation.Description readTranslationRotation(DataInput in)
	throws IOException {
		return new TranslationRotation.Description(
			in.readDouble(),
			in.readDouble(),
			new Vector3d(
				in.readDouble(),
				in.readDouble(),
				in.readDouble()
			)
		);
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
	 *
	 * Ideally, we'd use named static factory methods here, except we have a ton of final fields,
	 * and it would be a huge PITA to pass them all as constructor arguments,
	 * so this constructor signature hack will have to do for now.
	 */

	private int posPairIndex(int posi1, int posi2) {

		if (posi2 > posi1) {
			throw new IllegalArgumentException("posi2 should be strictly less than posi1");
		}

		return posi1*(posi1 - 1)/2 + posi2;
	}

	public IndicesStatic indicesStatic(int ffi) {
		return indicesStatic[ffi];
	}

	public IndicesSingle indicesSingles(int ffi, int posi, int confi) {

		// convert the conf index to a frag index
		int fragi = positions[posi].confs[confi].fragIndex;

		return indicesSinglesByFrag(ffi, posi, fragi);
	}

	public IndicesSingle indicesSinglesByFrag(int ffi, int posi, int fragi) {
		return indicesSingles[ffi][posi][fragi];
	}

	public IndicesPair indicesPairs(int ffi, int posi1, int confi1, int posi2, int confi2) {

		// convert the conf indices to frag indices
		int fragi1 = positions[posi1].confs[confi1].fragIndex;
		int fragi2 = positions[posi2].confs[confi2].fragIndex;

		return indicesPairsByFrags(ffi, posi1, fragi1, posi2, fragi2);
	}

	public IndicesPair indicesPairsByFrags(int ffi, int posi1, int fragi1, int posi2, int fragi2) {
		return indicesPairs[ffi][posPairIndex(posi1, posi2)][fragi1][fragi2];
	}

	public int getStaticAtomIndex(int atomi) {
		// easy peasy
		return atomi;
	}

	public int getConfAtomIndex(int posi, int atomi) {

		// just in case...
		assert (atomi < positions[posi].maxNumAtoms) :
			String.format("atom %d can't belong to position %d, which only has %d atoms", atomi, posi, positions[posi].maxNumAtoms);

		return confAtomOffsetsByPos[posi] + atomi;
	}

	public Integer findPosIndex(int atomi) {

		// look in the static range
		if (atomi < confAtomOffsetsByPos[0]) {
			return PosInter.StaticPos;
		}

		// look in the conf ranges
		for (int posi=0; posi<positions.length - 1; posi++) {
			if (atomi < confAtomOffsetsByPos[posi + 1]) {
				return posi;
			}
		}
		if (atomi < maxNumConfAtoms) {
			return positions.length - 1;
		}

		// not in any range
		return null;
	}

	/**
	 * Returns the average the number of atoms pairs for the given position interactions, over all conformations.
	 * Useful for predicting the cost of energy minimizations.
	 */
	public long avgAtomPairs(List<PosInter> inters) {
		long count = 0;
		for (var inter : inters) {
			for (var ffi=0; ffi<forcefieldIds.length; ffi++) {
				count += avgAtomPairs(ffi, inter);
			}
		}
		return count;
	}

	public long avgAtomPairs(int ffi, PosInter inter) {

		long count = 0;
		int n = 0;

		if (inter.posi1 == inter.posi2) {
			if (inter.posi1 == PosInter.StaticPos) {
				count += indicesStatic(ffi).size();
				n += 1;
			} else {
				for (int confi=0; confi<numConf(inter.posi1); confi++) {
					count += indicesSingles(ffi, inter.posi1, confi).sizeInternals();
					n += 1;
				}
			}
		} else if (inter.posi1 == PosInter.StaticPos) {
			for (int confi=0; confi<numConf(inter.posi2); confi++) {
				count += indicesSingles(ffi, inter.posi2, confi).sizeStatics();
				n += 1;
			}
		} else if (inter.posi2 == PosInter.StaticPos) {
			for (int confi=0; confi<numConf(inter.posi1); confi++) {
				count += indicesSingles(ffi, inter.posi1, confi).sizeStatics();
				n += 1;
			}
		} else {
			for (int confi1=0; confi1<numConf(inter.posi1); confi1++) {
				for (int confi2=0; confi2<numConf(inter.posi2); confi2++) {
					count += indicesPairs(ffi, inter.posi1, confi1, inter.posi2, confi2).size();
					n += 1;
				}
			}
		}

		// take the average
		return count/n;
	}

	public double[] ffparams(int ffi, int paramsi) {
		return ffparams[ffi][paramsi];
	}

	@Override
	public int countSingles() {
		int count = 0;
		for (int posi1=0; posi1<positions.length; posi1++) {
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

	public int numFrag(int posi) {
		return positions[posi].numFrags;
	}

	@Override
	public String name(int posi) {
		return positions[posi].name;
	}

	@Override
	public String confId(int posi, int confi) {
		return positions[posi].confs[confi].id;
	}

	@Override
	public String confType(int posi, int confi) {
		return positions[posi].confs[confi].type;
	}

	@Override
	public SeqSpace seqSpace() {
		return seqSpace;
	}

	@Override
	public String wildType(int posi) {
		return positions[posi].wildType;
	}

	@Override
	public boolean hasMutations(int posi) {
		return positions[posi].hasMutations;
	}

	/**
	 * Finds the first position with the given name, or null if no positions match.
	 * WARNING: position names are not guaranteed to be unique.
	 */
	public Pos findPos(String name) {
		for (Pos pos : positions) {
			if (pos.name.equals(name)) {
				return pos;
			}
		}
		return null;
	}

	/**
	 * Finds the first static atom with the given name and returns its index,
	 * or null if no static atom names match.
	 * WARNING: static atom names are not guaranteed to be unique.
	 */
	public Integer findStaticAtomIndex(String name) {
		for (int i=0; i<staticNames.length; i++) {
			if (staticNames[i].equals(name)) {
				return i;
			}
		}
		return null;
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
		AssignedCoords coords = new AssignedCoords(this, assignments);
		coords.copyCoords();
		coords.makeDofs();
		return coords;
	}

	@Override
	public int hashCode() {
		return hash;
	}
}
