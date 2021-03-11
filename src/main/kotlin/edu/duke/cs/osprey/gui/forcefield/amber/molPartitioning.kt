package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.*
import kotlin.NoSuchElementException
import kotlin.collections.ArrayList


// these residue types match the residue types used in AmberTools19 pdb4amber

private val aminoAcidResTypes = setOf(
	"ALA", "A",
	"ARG", "R",
	"ASN", "N",
	"ASP", "D", "ASH", "AS4",
	"CYS", "C", "CYM", "CYX",
	"GLU", "E", "GLH", "GL4",
	"GLN", "Q",
	"GLY", "G",
	"HIS", "H", "HIP", "HIE", "HID",
	"HYP",
	"ILE", "I",
	"LEU", "L",
	"LYS", "K", "LYN",
	"MET", "M",
	"PHE", "F",
	"PRO", "P",
	"SER", "S",
	"THR", "T",
	"TRP", "W",
	"TYR", "Y",
	"VAL", "V"
)

private val dnaResTypes = setOf(
	"DG", "GUA", "DG5", "DG3", "DGN",
	"DC", "CYT", "DC5", "DC3", "DCN", "DCP",
	"DA", "ADE", "DA5", "DA3", "DAN", "DAP",
	"DT", "THY", "DT5", "DT3"
)

private val rnaResTypes = setOf(
	"G", "GUA", "G5", "G3", "GN", "RG", "RG3", "RG5", "RGN", "GF2", "M2G", "YYG", "7MG", "OMG", "2MG",
	"C", "CYT", "CP", "C5", "C3", "CN", "RC", "RC5", "RC3", "RCN", "CFZ", "5MC", "OMC",
	"A", "ADE", "AP", "A5", "A3", "AN", "RA", "RA3", "RA5", "AF2", "1MA",
	"U", "URA", "U3", "U5", "UN", "RU", "RU3", "RU5", "RUN", "UFT", "5MU", "H2U", "PSU",
	"T", "THY", "T3", "T5", "TN", "RT", "RT3", "RT5", "RTN"
)

private val solventResTypes = setOf(
	"WAT", "HOH", "TIP3", "TIP4", "TIP5", "SPCE", "SPC", "SOL"
)

private val atomicIonResTypes = setOf(
	"Na+", "Li+", "Mg+", "Rb+", "MG", "Cs+", "POT", "SOD", "MG2",
	"CAL", "RUB", "LIT", "ZN2", "CD2", "NA", "K+", "K", "NA+",
	"Cl-", "Br-", "F-", "I-", "CLA", "CL", "BR", "CL-"
)

private val syntheticResTypes = setOf(
	// amber calls these "extra points"
	// also, they're album types
	"EP", "LP"
)

data class ForcefieldName(
	val name: String,
	val atomTypes: AmberAtomTypes? = null
) {

	val atomTypesOrThrow get() =
		atomTypes ?: throw NoSuchElementException("forcefield $name doesn't have atom types for Antechamber")

	override fun toString() = name

	companion object {

		val ff96 = ForcefieldName("ff96", AmberAtomTypes.Amber)
		val ff14SB = ForcefieldName("protein.ff14SB", AmberAtomTypes.Amber)

		val DNAOL15 = ForcefieldName("DNA.OL15", AmberAtomTypes.Amber)
		val RNAOL15 = ForcefieldName("RNA.OL3", AmberAtomTypes.Amber)
		val tip3p = ForcefieldName("water.tip3p", AmberAtomTypes.Amber)

		val gaff = ForcefieldName("gaff", AmberAtomTypes.Gaff)
		val gaff2 = ForcefieldName("gaff2", AmberAtomTypes.Gaff2)
	}
}

enum class MoleculeType(
	val isPolymer: Boolean,
	val forcefieldNames: List<ForcefieldName>
) {

	Protein(true, listOf(
		ForcefieldName.ff96, // dlab's favorite and time-tested protein forecfield
		ForcefieldName.ff14SB // currently recommended by AmberTools19
	)),

	DNA(true, listOf(
		ForcefieldName.DNAOL15 // currently recommended by AmberTools19
	)),

	RNA(true, listOf(
		ForcefieldName.RNAOL15 // currently recommended by AmberTools19
	)),

	Solvent(false, listOf(
		ForcefieldName.tip3p // currently recommended by AmberTools19
	)),

	AtomicIon(false, listOf(
		ForcefieldName.tip3p // currently recommended by AmberTools19
	)),

	Synthetic(false, listOf(
		// not real molecules, no forcefield needed
	)),

	SmallMolecule(false, listOf(
		ForcefieldName.gaff2, // currently recommended by AmberTools19
		ForcefieldName.gaff
	));


	// the default option is the first in the list
	val defaultForcefieldName get() = forcefieldNames.firstOrNull()

	val defaultForcefieldNameOrThrow get() =
		defaultForcefieldName ?: throw NoSuchElementException("molecule type $this has no default forcefield")

	companion object {
		operator fun get(resType: String) =
			when {

				// all the usual stuff
				aminoAcidResTypes.contains(resType) -> Protein
				dnaResTypes.contains(resType) -> DNA
				rnaResTypes.contains(resType) -> RNA
				solventResTypes.contains(resType) -> Solvent
				atomicIonResTypes.contains(resType) -> AtomicIon

				// we can typically ignore these, since they don't represent real molecules
				syntheticResTypes.contains(resType) -> Synthetic

				// assume anything we don't know is a small molecule
				else -> SmallMolecule
			}
	}
}

/**
 * Returns the types of the molecule
 * based on AMBER rules for residue classification.
 */
fun Molecule.findTypes(): Set<MoleculeType> {

	// for non-polymers, assume the whole molecule is a small molecule
	if (this !is Polymer) {
		return setOf(MoleculeType.SmallMolecule)
	}

	return chains
		.flatMap { it.residues }
		.map { MoleculeType[it.type] }
		.toSet()
}

/**
 * Returns the type of the molecule
 * based on AMBER rules for residue classification,
 * or throws an exception.
 */
fun Molecule.findTypeOrThrow(): MoleculeType {
	val types = findTypes()
	if (types.size == 1) {
		return types.first()
	}
	throw NoSuchElementException("No unique molecule type found in $types")
}

/**
 * Partition a single Molecule into a list of Molecules
 * based on AMBER rules for residue classification.
 *
 * Also returns a map between atoms in the input molecule (A side)
 * and atoms in the partitioned molecules (B side).
 */
fun Molecule.partitionAndAtomMap(combineSolvent: Boolean = true): Pair<List<Pair<MoleculeType,Molecule>>,AtomMap> {

	// for non-polymers, assume the whole molecule is a small molecule
	if (this !is Polymer) {
		return listOf(MoleculeType.SmallMolecule to this) to AtomMap.identity(atoms)
	}

	data class Partitioned(
		val moltype: MoleculeType,
		val chainId: String,
		val residues: List<Polymer.Residue>
	)

	// create the partition
	val partition = ArrayList<Partitioned>()
	for (chain in chains) {

		// chains in PDB files typically have lots of different molecules in them,
		// so separate them out by contiguous moltypes in PDB residue order
		var currentMoltype: MoleculeType? = null
		var residues: MutableList<Polymer.Residue>? = null
		for (res in chain.residues) {

			// if the moltype changed
			val moltype = MoleculeType[res.type]
			if (moltype != currentMoltype) {
				currentMoltype = moltype

				// start a new group of residues
				residues = ArrayList()
				residues.add(res)
				partition.add(Partitioned(moltype, chain.id, residues))

			} else {

				// otherwise, append to the existing group
				residues!!.add(res)
			}
		}
	}

	val combinedAtomMap = AtomMap()

	// create a molecule for each item in the partition
	var mols = partition.flatMap { (moltype, chainId, residues) ->

		if (moltype.isPolymer) {

			val atomMap = AtomMap()

			// map all the residues to a new polymer
			val polymer = Polymer("Chain $chainId")
			val chain = Polymer.Chain(chainId).also { polymer.chains.add(it) }

			// copy the atoms
			val srcAtoms = residues.flatMap { it.atoms }
			for (srcAtom in srcAtoms) {
				val dstAtom = srcAtom.copy()
				atomMap.add(srcAtom, dstAtom)
				polymer.atoms.add(dstAtom)
			}

			// copy the bonds (within the same partition item)
			for (srcAtom in srcAtoms) {
				val dstAtom = atomMap.getBOrThrow(srcAtom)
				for (srcBonded in bonds.bondedAtoms(srcAtom)) {
					val dstBonded = atomMap.getB(srcBonded)
						?: throw IllegalArgumentException(
							"$srcBonded is not in the residues in this $moltype polymer, but is bonded to $srcAtom which is"
						)
					polymer.bonds.add(dstAtom, dstBonded)
				}
			}

			// copy the residues
			for (res in residues) {
				chain.residues.add(Polymer.Residue(res.id, res.type, res.atoms.mapNotNull { atomMap.getB(it) }))
			}

			combinedAtomMap.addAll(atomMap)

			return@flatMap listOf(moltype to polymer)

		} else {

			// map each residue to a new molecule
			return@flatMap residues.map { res ->

				val mol = Molecule(res.type, res.type)

				val atomMap = AtomMap()

				// copy the atoms
				for (srcAtom in res.atoms) {
					val dstAtom = srcAtom.copy()
					atomMap.add(srcAtom, dstAtom)
					mol.atoms.add(dstAtom)
				}

				// copy the bonds (within the same residue)
				for (srcAtom in res.atoms) {
					val dstAtom = atomMap.getBOrThrow(srcAtom)
					for (srcBonded in bonds.bondedAtoms(srcAtom)) {
						val dstBonded = atomMap.getB(srcBonded) ?: throw IllegalArgumentException("bond spans multiple molecules")
						mol.bonds.add(dstAtom, dstBonded)
					}
				}
				// TODO: will that work for multi-residue non-protein/nucleic acid polymers? eg glycans?

				combinedAtomMap.addAll(atomMap)

				return@map moltype to mol
			}
		}
	}

	if (combineSolvent) {

		mols
			.filter { (type, _) -> type == MoleculeType.Solvent }
			.map { (_, mol) -> mol }
			.takeIf { it.isNotEmpty() }
			?.let { solventMols ->

				// combine all the solvent molecules into a "molecule"
				val name = solventMols
					.map { it.name }
					.toSet()
					.sorted()
					.joinToString(", ")
				val combinedSolvent = Molecule(name, solventMols.first().type)
				solventMols.forEachIndexed { i, solventMol ->
					combinedSolvent.atoms.addAll(solventMol.atoms)
				}

				mols = mols
					.filter { (type, _) -> type != MoleculeType.Solvent }
					.toMutableList() + listOf(MoleculeType.Solvent to combinedSolvent)
			}
	}

	return mols to combinedAtomMap
}

fun Molecule.partition(combineSolvent: Boolean = true): List<Pair<MoleculeType,Molecule>> =
	partitionAndAtomMap(combineSolvent).first

fun Collection<Molecule>.partition(combineSolvent: Boolean = true) =
	flatMap { it.partition(combineSolvent) }

fun Collection<Molecule>.partitionAndAtomMap(combineSolvent: Boolean = true): Pair<List<Pair<MoleculeType,Molecule>>,AtomMap> {
	val combinedMols = ArrayList<Pair<MoleculeType,Molecule>>()
	val combinedAtomMap = AtomMap()
	for (srcMol in this) {
		val (dstMol, atomMap) = srcMol.partitionAndAtomMap(combineSolvent)
		combinedMols.addAll(dstMol)
		combinedAtomMap.addAll(atomMap)
	}
	return combinedMols to combinedAtomMap
}