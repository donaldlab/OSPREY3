package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.molscope.tools.associateIdentity
import edu.duke.cs.osprey.gui.forcefield.AtomIndex
import edu.duke.cs.osprey.gui.io.*
import edu.duke.cs.osprey.gui.tools.CombineCollisionException
import edu.duke.cs.osprey.gui.tools.IntPair
import edu.duke.cs.osprey.gui.tools.combineMaps
import edu.duke.cs.osprey.service.services.ForcefieldParamsRequest
import edu.duke.cs.osprey.service.services.MoleculeFFInfoRequest
import edu.duke.cs.osprey.service.services.TypesRequest
import java.util.*


enum class AmberAtomTypes(val id: String) {
	Gaff("gaff"),
	Gaff2("gaff2"),
	Amber("amber"),
	BCC("bcc"),
	SYBYL("sybyl")
}

data class AmberTypes(
	val ffname: ForcefieldName,
	val atomTypes: Map<Int,String>,
	val bondTypes: Map<IntPair,String>,
	val atomCharges: Map<Int,String>
) {

	constructor(ffname: ForcefieldName, mol2Metadata: Mol2Metadata, atomIndex: AtomIndex) : this(
		ffname,
		mol2Metadata.atomTypes.mapKeys { (atom, _) ->
			atomIndex.getOrThrow(atom)
		},
		mol2Metadata.bondTypes.mapKeys { (bond, _) ->
			IntPair(atomIndex.getOrThrow(bond.a), atomIndex.getOrThrow(bond.b))
		},
		mol2Metadata.atomCharges.mapKeys { (atom, _) ->
			atomIndex.getOrThrow(atom)
		}
	)

	fun toMol2Metadata(mol: Molecule, atomIndex: AtomIndex): Mol2Metadata {

		val metadata = Mol2Metadata()

		// copy over the atom and bond types
		for ((atomi, type) in atomTypes) {
			metadata.atomTypes[atomIndex.getOrThrow(atomi)] = type
		}
		for ((bond, type) in bondTypes) {
			metadata.bondTypes[AtomPair(atomIndex.getOrThrow(bond.i1), atomIndex.getOrThrow(bond.i2))] = type
		}
		for ((atomi, charge) in atomCharges) {
			metadata.atomCharges[atomIndex.getOrThrow(atomi)] = charge
		}

		if (mol is Polymer) {
			for (chain in mol.chains) {
				for (res in chain.residues) {
					metadata.dictionaryTypes[res] = Mol2Metadata.defaultDictionaryType
				}
			}
		}

		return metadata
	}

	companion object {

		fun combine(
			types1: AmberTypes, preferredAtomIndices1: Set<Int>, ignoredAtomIndices1: Set<Int>,
			types2: AmberTypes, preferredAtomIndices2: Set<Int>, ignoredAtomIndices2: Set<Int>
		): AmberTypes {

			val ffname = if (types1.ffname == types2.ffname) {
				types1.ffname
			} else {
				throw IllegalArgumentException("""
					|two Amber types disagree on the forcefield name:
					|   ${types1.ffname}
					|   ${types2.ffname}
				""".trimMargin())
			}

			val combinedAtomTypes = try {
				combineMaps(
					types1.atomTypes, preferredAtomIndices1, ignoredAtomIndices1,
					types2.atomTypes, preferredAtomIndices2, ignoredAtomIndices2
				)
			} catch (ex: CombineCollisionException) {
				throw IllegalArgumentException("""
					|two Amber types disagree on types for atom ${ex.key}:
					|	${ex.val1}
					|	${ex.val2}
				""".trimMargin())
			}

			val combinedBondTypes = try {

				// prefer and ignore bonds types based on the preferred and ignored atoms
				fun Set<IntPair>.preferred(preferredAtomIndices: Set<Int>) =
					filter { (i1, i2) -> i1 in preferredAtomIndices && i2 in preferredAtomIndices }
					.toSet()
				fun Set<IntPair>.ignored(ignoredAtomIndices: Set<Int>) =
					filter { (i1, i2) -> i1 in ignoredAtomIndices || i2 in ignoredAtomIndices }
					.toSet()

				combineMaps(
					types1.bondTypes, types1.bondTypes.keys.preferred(preferredAtomIndices1), types1.bondTypes.keys.ignored(ignoredAtomIndices1),
					types2.bondTypes, types2.bondTypes.keys.preferred(preferredAtomIndices2), types2.bondTypes.keys.ignored(ignoredAtomIndices2)
				)
			} catch (ex: CombineCollisionException) {
				throw IllegalArgumentException("""
					|two Amber types disagree on types for bond ${ex.key}:
					|	${ex.val1}
					|	${ex.val2}
				""".trimMargin())
			}

			val combinedAtomCharges = try {
				combineMaps(
					types1.atomCharges, preferredAtomIndices1, ignoredAtomIndices1,
					types2.atomCharges, preferredAtomIndices2, ignoredAtomIndices2
				)
			} catch (ex: CombineCollisionException) {
				throw IllegalArgumentException("""
					|two Amber types disagree on charges for atom ${ex.key}:
					|	${ex.val1}
					|	${ex.val2}
				""".trimMargin())
			}

			return AmberTypes(ffname, combinedAtomTypes, combinedBondTypes, combinedAtomCharges)
		}
	}
}

enum class AmberChargeMethod(val id: String) {
	RESP("resp"),
	AM1BCC("bcc"),
	CM1("cm1"),
	CM2("cm2"),
	ESP("esp"),
	Mulliken("mul"),
	Gasteiger("gas")
}


data class AmberChargeGeneration(
	val method: AmberChargeMethod,
	val netCharge: Int,
	val minimizationSteps: Int
)

/**
 * Calculate Amber types for atoms and bonds.
 * Pass values to `chargeMethod` and `netCharge` to calculate partial charges for small molecules as well.
 */
suspend fun Molecule.calcTypesAmber(
	molType: MoleculeType,
	atomIndex: AtomIndex,
	ffname: ForcefieldName = molType.defaultForcefieldNameOrThrow,
	generateCharges: AmberChargeGeneration? = null
): AmberTypes {

	val dst = this

	val (srcMetadata, atomMap) = when (molType) {

		MoleculeType.SmallMolecule -> {

			// call osprey service with small molecule settings
			val request = TypesRequest.SmallMoleculeSettings(
				mol2 = dst.toMol2(),
				atomTypes = ffname.atomTypesOrThrow.id,
				chargeSettings = generateCharges?.let {
					TypesRequest.ChargeSettings(
						chargeMethod = it.method.id,
						netCharge = it.netCharge,
						numMinimizationSteps = it.minimizationSteps
					)
				}
			).toRequest()
			val (src, srcMetadata) = Molecule.fromMol2WithMetadata(OspreyService.types(request).mol2)

			// Tragically, we antechamber doesn't write the residue info back into the mol2 file,
			// so we can't use our usual MoleculeMapper to do the atom mapping here.
			// Thankfully, the atom order is preserved, so we can use that to do the mapping instead.
			val atomMap = src.atoms.zip(dst.atoms)
				.associateIdentity { (srcAtom, dstAtom) -> srcAtom to dstAtom }
			srcMetadata to atomMap
		}

		// for everything else, call osprey service with molecule settings
		else -> {

			val request = TypesRequest.MoleculeSettings(
				pdb = dst.toPDB(
					// amber will only see the disulfide bonds if we send
					// CYX residues, SSBOND records, and CONECT records
					translateSSasCYX = true,
					includeSSBondConect = true,
					// amber also needs histidine protonation state explicitly described
					translateHIStoEDP = true,
					// these errors will cause downstream problems, so fail loudly and early
					throwOnNonChainPolymerAtoms = true
				),
				ffname = ffname.name
			).toRequest()
			val (src, srcMetadata) = Molecule.fromMol2WithMetadata(OspreyService.types(request).mol2)

			// check for unmapped atoms with the dst->src mapping
			MoleculeMapper(dst, src)
				.let { dstToSrc ->
					dst.atoms.filter { dstToSrc.mapAtom(it) == null }
				}
				.takeIf { it.isNotEmpty() }
				?.map { atom ->
					// get descriptive info for the atom
					(dst as? Polymer)
						?.findChainAndResidue(atom)
						?.let { (chain, res) ->
							"${atom.name} @ ${chain.id}${res.id}"
						}
						?: atom.name
				}
				?.let { unmappedAtoms ->
					throw NoSuchElementException("LEaP didn't generate atom info for ${unmappedAtoms.size} atom(s):\n$unmappedAtoms")
				}

			// use the molecule mapper to map the metadata back with the src->dst mapping
			val srcToDst = MoleculeMapper(src, dst)
			srcMetadata to src.atoms.associateIdentity { srcAtom ->
				srcAtom to srcToDst.mapAtom(srcAtom)
			}
		}
	}

	// translate the metadata back
	// ignore atoms that aren't in the input molecules (sometimes amber adds stuff)
	// and check to make sure we have indices for all the atoms
	val dstMetadata = Mol2Metadata()
	for ((srcAtom, type) in srcMetadata.atomTypes) {
		val dstAtom = atomMap.getValue(srcAtom) ?: continue
		atomIndex.getOrThrow(dstAtom)
		dstMetadata.atomTypes[dstAtom] = type
	}
	for ((srcBond, type) in srcMetadata.bondTypes) {
		val dsta = atomMap.getValue(srcBond.a) ?: continue
		val dstb = atomMap.getValue(srcBond.b) ?: continue
		atomIndex.getOrThrow(dsta)
		atomIndex.getOrThrow(dstb)
		dstMetadata.bondTypes[AtomPair(dsta, dstb)] = type
	}
	for ((srcAtom, charge) in srcMetadata.atomCharges) {
		val dstAtom = atomMap.getValue(srcAtom) ?: continue
		atomIndex.getOrThrow(dstAtom)
		dstMetadata.atomCharges[dstAtom] = charge
	}

	return AmberTypes(ffname, dstMetadata, atomIndex)
}


suspend fun Molecule.calcModsAmber(types: AmberTypes, atomIndex: AtomIndex): String? {
	val mol2 = toMol2(types.toMol2Metadata(this, atomIndex))
	return OspreyService.moleculeFFInfo(MoleculeFFInfoRequest(mol2, types.ffname.name)).ffinfo
}

data class AmberMolParams(
	val mol: Molecule,
	val atomIndex: AtomIndex,
	val types: AmberTypes,
	val frcmods: List<String>
)

data class AmberParams(
	val top: String,
	val crd: String
)

suspend fun List<AmberMolParams>.calcParamsAmber(): AmberParams {

	val response = OspreyService.forcefieldParams(ForcefieldParamsRequest(
		map {
			ForcefieldParamsRequest.MolInfo(
				mol2 = it.mol.toMol2(it.types.toMol2Metadata(it.mol, it.atomIndex)),
				ffname = it.types.ffname.name,
				ffinfos = it.frcmods
			)
		}
	))

	return AmberParams(response.params, response.coords)
}

suspend fun AmberMolParams.calcParamsAmber() =
	listOf(this).calcParamsAmber()
