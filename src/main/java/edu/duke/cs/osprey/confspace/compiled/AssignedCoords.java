package edu.duke.cs.osprey.confspace.compiled;


import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import org.joml.Vector3d;

import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Supplier;


/**
 * A copy of the conf space atom coords with the desired assignments.
 */
public class AssignedCoords {

	public final ConfSpace confSpace;

	/** conf indices for the design positions, in order */
	public final int[] assignments;

	/** atom coords for the static atoms and the conformation atoms */
	public final CoordsList coords;

	/** degrees of freedom that modify the atom coords */
	public final List<DegreeOfFreedom> dofs = new ArrayList<>();

	public AssignedCoords(ConfSpace confSpace, int[] assignments) {

		this.confSpace = confSpace;
		this.assignments = assignments;

		coords = new CoordsList(confSpace.maxNumConfAtoms);
	}

	public void copyCoords() {

		// copy over the static atoms first, so DoFs can modify them
		coords.copyFrom(confSpace.staticCoords, 0);

		// then copy the conformation coords
		for (ConfSpace.Pos pos : confSpace.positions) {

			// get the conf, or skip this position if nothing was assigned
			int confi = assignments[pos.index];
			if (confi == ConfSpace.NotAssigned) {
				continue;
			}
			ConfSpace.Conf conf = pos.confs[confi];

			// copy over the coords
			coords.copyFrom(conf.coords, confSpace.confAtomOffsetsByPos[pos.index]);
		}
	}

	public void makeDofs() {

		dofs.clear();

		// first, make the molecule motions and convert them into degrees of freedom
		for (int moli=0; moli<confSpace.molInfos.length; moli++) {
			ConfSpace.MolInfo molInfo = confSpace.molInfos[moli];
			for (int i=0; i<molInfo.motions.length; i++) {
				// TODO: do we always need this motion?
				ContinuousMotion motion = molInfo.motions[i].build(this, moli);
				motion.appendDofs(dofs);
			}
		}

		// then, make motions for conformations
		for (ConfSpace.Pos pos : confSpace.positions) {

			// get the conf, or skip this position if nothing was assigned
			int confi = assignments[pos.index];
			if (confi == ConfSpace.NotAssigned) {
				continue;
			}
			ConfSpace.Conf conf = pos.confs[confi];

			// make the motions and convert them into degrees of freedom
			for (int i=0; i<conf.motions.length; i++) {
				ContinuousMotion motion = conf.motions[i].build(this, pos);
				motion.appendDofs(dofs);
			}
		}

		// just in case...
		assert (dofs.size() <= confSpace.maxNumDofs);
	}

	public double getStaticEnergy(int ffi) {
		return confSpace.staticEnergies[ffi];
	}

	public double getInternalEnergy(int ffi, int posi) {

		// get the assignment, or 0 if nothing was assigned
		int confi = assignments[posi];
		if (confi == ConfSpace.NotAssigned) {
			return 0.0;
		}

		return confSpace.positions[posi].confs[confi].energies[ffi];
	}

	public ConfSpace.IndicesStatic getIndices(int ffi) {
		return confSpace.indicesStatic(ffi);
	}

	public ConfSpace.IndicesSingle getIndices(int ffi, int posi) {

		// get the assignment, or null if nothing was assigned
		int confi = assignments[posi];
		if (confi == ConfSpace.NotAssigned) {
			return null;
		}

		return confSpace.indicesSingles(ffi, posi, confi);
	}

	public ConfSpace.IndicesPair getIndices(int ffi, int posi1, int posi2) {

		// get the assignments, or null if nothing was assigned
		int confi1 = assignments[posi1];
		int confi2 = assignments[posi2];
		if (confi1 == ConfSpace.NotAssigned || confi2 == ConfSpace.NotAssigned) {
			return null;
		}

		return confSpace.indicesPairs(ffi, posi1, confi1, posi2, confi2);
	}

	public double[] getParams(int ffi, int paramsi) {
		return confSpace.ffparams(ffi, paramsi);
	}

	public int getStaticIndex(int atomi) {
		return confSpace.getStaticAtomIndex(atomi);
	}

	public int getConfIndex(int posi, int atomi) {
		return confSpace.getConfAtomIndex(posi, atomi);
	}

	/**
	 * Converts the coordinates to a Molecule instance.
	 * The molecule will have no bonds though,
	 * since bond information is not available here.
	 */
	public Molecule toMol() {

		// NOTE: This function doesn't get called too much,
		// so it doesn't have to be super fast.
		// It just needs to work correctly.

		Function<Integer,String> makeResKey = resInfoIndex -> {
			var resInfo = confSpace.resInfos[resInfoIndex];
			return String.format("%s%s", resInfo.chainId, resInfo.resId);
		};

		// figure out which res infos to use with the current assignments
		// start with res infos from the static atoms
		Map<String,Integer> staticResInfos = new HashMap<>();
		for (int i=0; i<confSpace.numStaticAtoms; i++) {
			var resInfoIndex = confSpace.staticResInfoIndices[i];
			staticResInfos.put(makeResKey.apply(resInfoIndex), resInfoIndex);
		}

		// then get the res infos from conf atoms
		Map<String,Integer> assignedResInfos = new HashMap<>();
		for (int posi=0; posi<confSpace.numPos(); posi++) {

			// get the conf
			int confi = assignments[posi];
			if (confi == Conf.Unassigned) {
				continue;
			}
			ConfSpace.Conf conf = confSpace.positions[posi].confs[confi];

			for (int atomi=0; atomi<conf.numAtoms; atomi++) {
				var resInfoIndex = conf.atomResInfoIndices[atomi];

				// check for collisions first
				// the compiled conf spaces allow making mutations that can't be expressed in PDB files
				// so converting to the old molecule class, which is modeled closely to the PDB format, causes problems
				var resKey = makeResKey.apply(resInfoIndex);
				var oldResInfoIndex = assignedResInfos.get(resKey);
				if (oldResInfoIndex != null && oldResInfoIndex != resInfoIndex) {
					throw new IllegalArgumentException(
						"Multiple assigments conflict at residue " + resKey + ", including position " + confSpace.positions[posi].name + "."
						+ " Can't pick unique residue type for residue-centric molecule formats, like PDB."
						+ " Probably caused by multiple design positions affecting the same residue."
					);
				}

				// no collisions, make the assignment
				assignedResInfos.put(resKey, resInfoIndex);
			}
		}

		// combine the static and assigned res info maps to make a single authoritative source for residue assignments
		Map<String,Integer> resInfoIndexMap = new HashMap<>();
		resInfoIndexMap.putAll(staticResInfos);
		resInfoIndexMap.putAll(assignedResInfos);

		// collect all the atoms by molecule
		// and also by chain and residue if possible
		class AtomInfo {
			final Vector3d pos = new Vector3d();
			public String name;
			int molInfoIndex;
			int resInfoIndex;
		}

		Map<Integer,Map<String,Map<Integer,List<AtomInfo>>>> residuedAtomInfos = new LinkedHashMap<>();
		Map<Integer,List<AtomInfo>> molAtomInfos = new LinkedHashMap<>();

		// adds atom infos to the right spot in the collections
		Consumer<AtomInfo> addInfo = info -> {
			if (info.resInfoIndex >= 0) {
				residuedAtomInfos
					.computeIfAbsent(info.molInfoIndex, index -> new LinkedHashMap<>())
					.computeIfAbsent(confSpace.resInfos[info.resInfoIndex].chainId, id -> new LinkedHashMap<>())
					.computeIfAbsent(info.resInfoIndex, index -> new ArrayList<>())
					.add(info);
			} else {
				molAtomInfos
					.computeIfAbsent(info.molInfoIndex, index -> new ArrayList<>())
					.add(info);
			}
		};

		// first, collect all the conf atoms
		for (int posi=0; posi<confSpace.numPos(); posi++) {

			// get the conf
			int confi = assignments[posi];
			if (confi == Conf.Unassigned) {
				continue;
			}
			ConfSpace.Conf conf = confSpace.positions[posi].confs[confi];

			for (int atomi=0; atomi<conf.numAtoms; atomi++) {

				AtomInfo atomInfo = new AtomInfo();
				coords.get(getConfIndex(posi, atomi), atomInfo.pos);
				atomInfo.name = conf.atomNames[atomi];
				atomInfo.molInfoIndex = conf.atomMolInfoIndices[atomi];
				atomInfo.resInfoIndex = conf.atomResInfoIndices[atomi];

				// just in case...
				if (atomInfo.resInfoIndex != resInfoIndexMap.get(makeResKey.apply(atomInfo.resInfoIndex))) {
					throw new Error("residue assignments from conf atoms don't match previously-determined residue assignments... from the conf atoms. Obviously this is a bug");
				}

				addInfo.accept(atomInfo);
			}
		}

		// collect all the static atoms
		for (int i=0; i<confSpace.numStaticAtoms; i++) {

			AtomInfo atomInfo = new AtomInfo();
			coords.get(getStaticIndex(i), atomInfo.pos);
			atomInfo.name = confSpace.staticNames[i];
			atomInfo.molInfoIndex = confSpace.staticMolInfoIndices[i];
			atomInfo.resInfoIndex = confSpace.staticResInfoIndices[i];

			// apply the residue assignments
			atomInfo.resInfoIndex = resInfoIndexMap.get(makeResKey.apply(atomInfo.resInfoIndex));

			addInfo.accept(atomInfo);
		}

		// finally, build the new molecule instance
		Molecule mol = new Molecule();

		Set<String> usedChainIds = new HashSet<>();
		Supplier<String> makeChainId = () -> {
			for (char id='A'; id<='Z'; id++) {
				String str = Character.toString(id);
				if (!usedChainIds.contains(str)) {
					usedChainIds.add(str);
					return str;
				}
			}
			throw new NoSuchElementException("can't find unused chain id");
		};

		Set<String> usedResIds = new HashSet<>();
		Supplier<String> makeResId = () -> {
			for (int id=1; id<9999; id++) {
				String str = Integer.toString(id);
				if (!usedResIds.contains(str)) {
					usedResIds.add(str);
					return str;
				}
			}
			throw new NoSuchElementException("can't find unused chain id");
		};

		BiConsumer<String,List<AtomInfo>> makeRes = (name, atomInfos) -> {

			// make the atoms
			ArrayList<Atom> atoms = new ArrayList<>();
			double[] coords = new double[atomInfos.size()*3];
			for (int i=0; i<atomInfos.size(); i++) {
				AtomInfo atomInfo = atomInfos.get(i);

				atoms.add(new Atom(atomInfo.name));

				coords[i*3] = atomInfo.pos.x;
				coords[i*3 + 1] = atomInfo.pos.y;
				coords[i*3 + 2] = atomInfo.pos.z;
			}

			mol.residues.add(new Residue(atoms, coords, name, mol));
		};

		// do the residued atoms first
		residuedAtomInfos.forEach((molInfoIndex, byChain) -> {
			byChain.forEach((chainId, byRes) -> {
				byRes.keySet().stream()
					.sorted(Comparator.comparing(resInfoIndex -> confSpace.resInfos[resInfoIndex].indexInChain))
					.forEach(resInfoIndex -> {

						List<AtomInfo> atomInfos = byRes.get(resInfoIndex);

						// make the residue
						ConfSpace.ResInfo resInfo = confSpace.resInfos[resInfoIndex];
						String name = String.format("%3s%2s%4s",
							resInfo.resType,
							resInfo.chainId,
							resInfo.resId
						);

						usedChainIds.add(resInfo.chainId);
						usedResIds.add(resInfo.resId);

						makeRes.accept(name, atomInfos);
					});
			});
		});

		// then do the molecule atoms
		molAtomInfos.forEach((molInfoIndex, atomInfos) -> {
			ConfSpace.MolInfo molInfo = confSpace.molInfos[molInfoIndex];

			// make one residue for all the atoms
			String name = String.format("%3s%2s%4s",
				molInfo.type != null ? molInfo.type : "XXX",
				makeChainId.get(),
				makeResId.get()
			);
			makeRes.accept(name, atomInfos);
		});

		return mol;
	}

	public String getAtomName(int atomi) {
		int posi = confSpace.findPosIndex(atomi);
		if (posi == PosInter.StaticPos) {
			return String.format("Static:%s",
				confSpace.staticNames[atomi]
			);
		} else {
			var pos = confSpace.positions[posi];
			int confi = assignments[posi];
			var conf = pos.confs[confi];
			int offset = atomi - confSpace.getConfAtomIndex(posi, 0);
			return String.format("%s:%s", pos.name, conf.atomNames[offset]);
		}
	}


	public static class AtomInfo {
		public final int posi;
		public final int atomi;
		public final int coordsi;
		public AtomInfo(int posi, int atomi, int coordsi) {
			this.posi = posi;
			this.atomi = atomi;
			this.coordsi = coordsi;
		}
	}

	/** returns information about the atom coords that are actually populated for this conformation */
	public List<AtomInfo> atoms() {
		var infos = new ArrayList<AtomInfo>();

		// add the static atoms
		for (int atomi=0; atomi<confSpace.staticCoords.size; atomi++) {
			infos.add(new AtomInfo(-1, atomi, confSpace.getStaticAtomIndex(atomi)));
		}

		// add the position atoms
		for (int posi=0; posi<confSpace.positions.length; posi++) {
			var pos = confSpace.positions[posi];
			var conf = pos.confs[assignments[posi]];
			for (int atomi=0; atomi<conf.numAtoms; atomi++) {
				infos.add(new AtomInfo(posi, atomi, confSpace.getConfAtomIndex(posi, atomi)));
			}
		}

		return infos;
	}
}
