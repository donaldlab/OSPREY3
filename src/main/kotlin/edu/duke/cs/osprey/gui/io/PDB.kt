package edu.duke.cs.osprey.gui.io

import cuchaz.kludge.tools.toString
import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.structure.Molecule as OspreyMolecule
import edu.duke.cs.osprey.structure.Atom as OspreyAtom
import edu.duke.cs.osprey.structure.Residue as OspreyResidue
import edu.duke.cs.osprey.structure.PDBIO
import edu.duke.cs.osprey.gui.prep.Proteins
import java.util.*
import kotlin.math.min


/**
 * Save the molecule to PDB format.
 */
fun Molecule.toPDB(
	/**
	 * If a Polymer contains atoms not in a chain, they will not be copied into the Osprey molecule.
	 * If true, throw an error when this happens.
	 * If false, print a warning to stderr.
	 *
	 * To avoid this issue, combine() your molecules with a ChainGenerator to generate chains for non-chain atoms.
	 */
	throwOnNonChainPolymerAtoms: Boolean = false,
	// TODO: once we're sure all the existing code doesn't have this issue,
	//  set this default to true to prevent future code from re-creating this issue
	translateSSasCYX: Boolean = false,
	comment: String? = null,
	energy: Double? = null,
	includeTer: Boolean = false,
	includeSSBondConect: Boolean = false,
	translateHIStoEDP: Boolean = false
): String {
	return PDBIO.write(
		toOspreyMol(throwOnNonChainPolymerAtoms, translateSSasCYX, translateHIStoEDP),
		comment,
		energy,
		includeTer,
		includeSSBondConect
	)
}


/**
 * Read a molecule in PDB format.
 */
fun Molecule.Companion.fromPDB(pdb: String): Molecule {
	return PDBIO.read(pdb).toMolecule("PDB")
}


/**
 * Convert an OSPREY molecule to a Molscope molecule.
 */
fun OspreyMolecule.toMolecule(name: String? = null): Molecule {

	val omol = this

	// if there's more than one residue, assume it's a polymer
	val mol = if (omol.residues.size > 1) {
		Polymer(name ?: omol.name)
	} else {
		Molecule(name ?: omol.name, omol.residues.firstOrNull()?.type)
	}

	// convert the atoms
	val atomMap = IdentityHashMap<OspreyAtom,Atom>()
	for (ores in omol.residues) {
		for (oatom in ores.atoms) {
			val atom = Atom(
				Element[oatom.elementNumber],
				oatom.name,
				oatom.coords[0],
				oatom.coords[1],
				oatom.coords[2]
			)
			atomMap[oatom] = atom
			mol.atoms.add(atom)
		}
	}

	// convert the bonds
	for (ores in omol.residues) {
		for (oatom in ores.atoms) {
			val atom = atomMap[oatom]!!
			for (bondedOatom in oatom.bonds) {
				val bondedAtom = atomMap[bondedOatom]!!
				mol.bonds.add(atom, bondedAtom)
			}
		}
	}

	if (mol is Polymer) {

		// find a unique chain id we can use
		val availableChainId = lazy {

			val chainIds = omol.residues
				.map { it.chainId }
				.toSet()
				.filter { !it.isWhitespace() }

			('A'..'Z')
				.firstOrNull { it !in chainIds }
				?.toString()
				?: throw NoSuchElementException("can't find a unique chain ID, IDs already in use: ${chainIds.sorted().joinToString(",")}")
		}

		// empty chain IDs should be replaced with something readable
		fun Char.translateChainId(): String =
			if (isWhitespace()) {
				availableChainId.value
			} else {
				toString()
			}

		// convert the residues
		for (res in omol.residues) {

			val ochainId = res.chainId.translateChainId()

			// get/make the chain
			val chain = mol.chains
				.find { it.id == ochainId }
				?: run {
					Polymer.Chain(ochainId).apply {
						mol.chains.add(this)
					}
				}

			chain.residues.add(Polymer.Residue(
				res.fullName.substring(5).trim(),
				res.type,
				res.atoms.map { atomMap[it]!! }
			))
		}
	}

	return mol
}


/**
 * Convert a Molscope molecule to an OSPREY molecule.
 */
fun Molecule.toOspreyMol(
	/**
	 * If a Polymer contains atoms not in a chain, they will not be copied into the Osprey molecule.
	 * If true, throw an error when this happens.
	 * If false, print a warning to stderr.
	 *
	 * To avoid this issue, combine() your molecules with a ChainGenerator to generate chains for non-chain atoms.
	 */
	throwOnNonChainPolymerAtoms: Boolean = false,
	// TODO: once we're sure all the existing code doesn't have this issue,
	//  set this default to true to prevent future code from re-creating this issue
	translateSSasCYX: Boolean = false,
	/**
	 * Amber can't read Histidine residues unless they're named based on their protonation state.
	 * If true, rename HIS residues to HIE, HID, or HIP depending on the protonation state
	 */
	translateHIStoEDP: Boolean = false
): OspreyMolecule {

	val mol = this
	val omol = OspreyMolecule()
	omol.name = mol.name

	val atomMap = IdentityHashMap<Atom,OspreyAtom>()

	fun Atom.toOsprey(): Pair<OspreyAtom,DoubleArray> {
		val oatom = OspreyAtom(name, element.symbol)
		atomMap[this] = oatom
		return oatom to doubleArrayOf(pos.x, pos.y, pos.z)
	}

	fun String.first(len: Int) = substring(0, min(len, length))

	// build the residue "full name", eg "ASN A  23"
	fun fullName(chainId: String, resId: String, resType: String) =
		"%3s%2s%4s".format(
			resType.first(3),
			chainId[0],
			resId.first(4)
		)

	when (mol) {
		is Polymer -> {

			val atomsCopied = Atom.identitySet()

			// put the atoms in multiple residues
			for (chain in mol.chains) {
				for (res in chain.residues) {

					val atoms = ArrayList<OspreyAtom>()
					val coords = ArrayList<DoubleArray>()
					for (atom in res.atoms) {
						val (oatom, oatomCoords) = atom.toOsprey()
						atoms.add(oatom)
						coords.add(oatomCoords)
					}

					// if desired, and if a disulfide bond is present, translate CYS residues to CYX residues
					val resType = if (translateSSasCYX && res.type == "CYS" && Proteins.isSSBonded(mol, res)) {
						"CYX"

					// if desired, translate HIS residues to describe the protonation state
					} else if (translateHIStoEDP && res.type == "HIS") {

						// determine the protonation state:
						// HIE has HE2
						// HID has HD1
						// HIP has both
						val hasHe2 = res.atoms.any { it.name.toLowerCase() == "he2" }
						val hasHd1 = res.atoms.any { it.name.toLowerCase() == "hd1" }
						if (hasHe2 && hasHd1) {
							"HIP"
						} else if (hasHe2) {
							"HIE"
						} else if (hasHd1) {
							"HID"
						} else {
							throw IllegalArgumentException("Unable to determine Histidine protonation state, no HE2 or HD1")
						}

					} else {
						res.type
					}

					omol.residues.add(OspreyResidue(atoms, coords, fullName(chain.id, res.id, resType), omol))

					atomsCopied.addAll(res.atoms)
				}
			}

			// warn about any non-chain atoms in the polymer
			val missingAtoms = atoms.size - atomsCopied.size
			if (missingAtoms > 0) {

				val msg = StringBuilder()
				msg.append("Polymer has $missingAtoms atoms that were not in chains")
				msg.append(", and not copied to the Osprey molecule.")

				// find the missing atoms
				val chainAtoms = Atom.identitySet()
				for (chain in mol.chains) {
					for (res in chain.residues) {
						chainAtoms.addAll(res.atoms)
					}
				}
				for (atom in mol.atoms) {
					if (atom !in chainAtoms) {
						msg.append("\n\t$atom @ ${atom.pos.toString(6)}")
					}
				}

				if (throwOnNonChainPolymerAtoms) {
					throw IllegalArgumentException(msg.toString())
				} else {
					System.err.println("WARNING: $msg")
				}
			}
		}
		else -> {

			// put all the atoms in a single residue
			val atoms = ArrayList<OspreyAtom>()
			val coords = ArrayList<DoubleArray>()
			for (atom in mol.atoms) {
				val (oatom, oatomCoords) = atom.toOsprey()
				atoms.add(oatom)
				coords.add(oatomCoords)
			}

			val resType = mol.type ?: mol.name.first(3).toUpperCase()
			omol.residues.add(OspreyResidue(atoms, coords, fullName("A", "1", resType), omol))
		}
	}

	// transfer the bonds
	for ((atom, oatom) in atomMap) {
		for (bonded in mol.bonds.bondedAtoms(atom)) {
			val obonded = atomMap[bonded]
				?: throw NoSuchElementException("bonded atom $bonded (bonded to $atom) not in this molecule: $mol")
			if (obonded !in oatom.bonds) {
				oatom.addBond(obonded)
			}
		}
	}

	return omol
}
