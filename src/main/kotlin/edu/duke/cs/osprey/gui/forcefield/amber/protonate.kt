package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.gui.io.OspreyService
import edu.duke.cs.osprey.gui.io.fromMol2
import edu.duke.cs.osprey.gui.io.toMol2
import edu.duke.cs.osprey.gui.io.toPDB
import edu.duke.cs.osprey.service.services.ProtonateRequest
import edu.duke.cs.osprey.service.services.ProtonationRequest
import kotlinx.coroutines.runBlocking


/**
 * Get all the supported protonation states for the specified atom.
 */
fun Molecule.protonations(atom: Atom): List<Protonation> {

	// get the number of bonded heavy atoms
	val numHeavy = bonds.bondedAtoms(atom)
		.filter { it.element != Element.Hydrogen }
		.size

	// get all the protonations supported for this atom
	return protonations.keys
		.filter { it.element == atom.element && it.numHeavy == numHeavy }
}

/**
 * Removes all bonded hydrogen atoms from the atom.
 */
fun Molecule.deprotonate(atom: Atom) {

	// remove all the hydrogens
	bonds.bondedAtoms(atom)
		.filter { it.element == Element.Hydrogen }
		.forEach {
			atoms.remove(it)
			if (this is Polymer) {
				findResidue(it)?.atoms?.remove(it)
			}
		}
}

/**
 * Removes all hydrogen atoms from the molecule.
 */
fun Molecule.deprotonate() {
	atoms
		.filter { it.element == Element.Hydrogen }
		.forEach {
			atoms.remove(it)
			if (this is Polymer) {
				findResidue(it)?.atoms?.remove(it)
			}
		}
}

fun Molecule.findProtonation(atom: Atom, numH: Int, hybridization: Hybridization): Protonation? =
	protonations(atom)
		.find { it.numH == numH && it.hybridization == hybridization }

fun Molecule.protonateBlocking(atom: Atom, protonation: Protonation) =
	runBlocking { protonate(atom, protonation) }


/**
 * Adds hydrogens atoms to the atom in the supplied protonation state.
 */
suspend fun Molecule.protonate(atom: Atom, protonation: Protonation) {

	val mol = this

	// get the GAFF atom types
	val amberTypes = protonations[protonation]
		?: throw IllegalArgumentException("protonation not supported for atom ${atom.name}: $protonation")

	// make a small molecule for just this atom and its heavy neighbors
	val smol = Molecule(mol.name)
		.apply {
			atoms.add(atom)
			mol.bonds.bondedAtoms(atom)
				.filter { it.element != Element.Hydrogen }
				.forEach {
					atoms.add(it)
					bonds.add(atom, it)
				}
		}
		.toMol2()

	// get the hydrogens whose names we shouldn't match
	val domainHydrogens = HashMap<String,Int?>().apply {
		if (mol is Polymer) {
			mol.findResidueOrThrow(atom).atoms
		} else {
			mol.atoms
		}
		.filter { it.element == Element.Hydrogen }
		.forEach { atom ->
			this[atom.name] = atom.name
				.filter { it.isDigit() }
				.takeIf { it.isNotBlank() }
				?.toInt()
		}
	}

	// pick a number for the hydrogen
	// try not to match other hydrogens in the molecule
	fun pickHydrogenName(i: Int): String {
		val atomNumber = atom.name.filter { it.isDigit() }
		var hNumber = "$atomNumber${i + 1}".toInt()
		if (hNumber in domainHydrogens.values) {
			hNumber = (domainHydrogens.values.maxBy { it ?: 0 } ?: 0) + 1
		}
		val hName = "H$hNumber"
		domainHydrogens[hName] = hNumber

		return hName
	}

	// build the service request

	val request = ProtonateRequest(
		mol2 = smol,
		atomName = atom.name,
		atomType = amberTypes.heavy,
		bonds = if (protonation.numHeavy == 1 && amberTypes.bondedHeavy != null && amberTypes.heavyBond != null) {

				// yup, try to set the bond and atoms types explicitly so we get the desired geometry
				// (leap doesn't recognize Sp hybridization explicitly,
				// and sometimes gets Sp2 hybridization wrong when we don't label the bonds correctly)

				// what's the other heavy atom?
				val heavyAtom = mol.bonds.bondedAtoms(atom)
					.filter { it.element != Element.Hydrogen }
					.takeIf { it.size == 1 }
					?.first()
					?: throw IllegalArgumentException("couldn't find unique heavy bonded atom")

				// set the other heavy atom type
				val heavyAtomType = amberTypes.bondedHeavy[heavyAtom.element]
					?: throw Error("protonation is Sp hybridized, but has no bonded heavy atom type for ${heavyAtom.element}")

				listOf(ProtonateRequest.Bond(
					atomName = heavyAtom.name,
					atomType = heavyAtomType,
					bondType = amberTypes.heavyBond
				))
			} else {
				emptyList()
			},
		hydrogens = (0 until protonation.numH)
			.map { i ->
				ProtonateRequest.Hydrogen(
					atomName = pickHydrogenName(i),
					atomType = amberTypes.hydrogen
				)
			}
	)

	val hmol = Molecule.fromMol2(OspreyService.protonate(request).mol2)

	// read the output mol2 and copy the new hydrogens
	val centerAtom = hmol.atoms.find { it.name == atom.name }
		?: throw Error("can't find central atom in LEaP output molecule")
	mol.deprotonate(atom)
	val res = (mol as? Polymer)?.findResidue(atom)
	hmol.bonds.bondedAtoms(centerAtom)
		.filter { it.element == Element.Hydrogen }
		.forEach {
			mol.atoms.add(it)
			mol.bonds.add(atom, it)
			res?.atoms?.add(it)
		}
}


data class ProtonatedAtom(
	val mol: Molecule,
	val chain: Polymer.Chain?,
	val res: Polymer.Residue?,
	val heavy: Atom,
	val light: Atom
) {

	fun included(): Boolean =
		mol.atoms.any { it === light }

	fun add() {

		// don't add more than once
		if (included()) {
			return
		}

		mol.atoms.add(light)
		res?.atoms?.add(light)
		mol.bonds.add(heavy, light)
	}

	fun remove() {
		mol.atoms.remove(light)
		res?.atoms?.remove(light)
		mol.bonds.remove(heavy, light)
	}

	/** a friendly description for the location of the atom */
	val location: String? get() =
		if (res != null && chain != null) {
			"${chain.id}${res.id}"
		} else {
			null
		}

	override fun toString() =
		location
			?.let { "${light.name} - ${heavy.name} @ $it" }
			?: "${light.name} - ${heavy.name}"
}


fun Molecule.inferProtonationBlocking() =
	runBlocking { inferProtonation() }

/**
 * Returns a list of heavy-hydrogen atom pairs based on inferred
 * forcefield atom and bond types.
 */
suspend fun Molecule.inferProtonation(): List<ProtonatedAtom> {

	val dst = this
	val dstAtoms = ArrayList<ProtonatedAtom>()

	// treat each molecule in the partition with the appropriate protocol
	val (partition, atomMap) = dst.partitionAndAtomMap(combineSolvent = true)
	partition@for ((type, src) in partition) {

		// TODO: allow user to pick the forcefields?
		val request = when (type) {

			// treat molecules with either leap or antechamber
			MoleculeType.Protein,
			MoleculeType.DNA,
			MoleculeType.RNA,
			MoleculeType.Solvent ->
				ProtonationRequest(src.toPDB(), type.defaultForcefieldNameOrThrow.name, null)

			MoleculeType.SmallMolecule -> {
				val ffname = type.defaultForcefieldNameOrThrow
				ProtonationRequest(src.toPDB(), ffname.name, ffname.atomTypesOrThrow.id)
			}

			// atomic ions don't have protonation
			MoleculeType.AtomicIon -> continue@partition

			// synthetics aren't real molecules, just ignore them
			MoleculeType.Synthetic -> continue@partition
		}
		val protonatedMol = Molecule.fromMol2(OspreyService.protonation(request).mol2)

		// translate the atoms to the input mol
		val srcAtoms = protonatedMol.translateHydrogens(src)
		for ((srcHeavy, h) in srcAtoms) {
			val dstHeavy = atomMap.getAOrThrow(srcHeavy)
			val chainRes = (dst as? Polymer)?.findChainAndResidue(dstHeavy)
			dstAtoms.add(ProtonatedAtom(
				mol = dst,
				chain = chainRes?.first,
				res = chainRes?.second,
				heavy = dstHeavy,
				light = h
			))
		}
	}

	return dstAtoms
}

private fun Molecule.translateHydrogens(dst: Molecule): List<Pair<Atom,Atom>> {

	val src: Molecule = this

	val mapper = src.mapTo(dst)

	// map the heavy,hydrogen pairs from src to dst
	return src.atoms
		.filter { it.element == Element.Hydrogen }
		.mapNotNull { srcH ->

			// get the bonded heavy atom
			val srcHeavy = src.bonds
				.bondedAtoms(srcH)
				.firstOrNull { it.element != Element.Hydrogen }
				?: return@mapNotNull null

			val dstHeavy = mapper.mapAtom(srcHeavy) ?: return@mapNotNull null

			return@mapNotNull dstHeavy to srcH
		}
}


data class Protonation(
	val element: Element,
	val numHeavy: Int,
	val numH: Int,
	val hybridization: Hybridization
)

private data class ProtonationTypes(
	val heavy: String,
	val hydrogen: String,
	val heavyBond: String? = null,
	val bondedHeavy: Map<Element,String>? = emptyMap()
)

private val protonations = mapOf(

	// NOTE: not sure if all these hybridizations are correct
	// they really only matter for Nitrogen though, when multiple protonations have the same number of hydrogen atoms
	// and the only way to distinguish them is with different hybridizations

	Protonation(Element.Carbon, 0, 2, Hybridization.Sp2) to ProtonationTypes("c1", "hc"), // eg methylene
	Protonation(Element.Carbon, 0, 3, Hybridization.Sp2) to ProtonationTypes("c2", "hc"), // eg methyl cation
	Protonation(Element.Carbon, 0, 4, Hybridization.Sp3) to ProtonationTypes("c3", "hc"), // eg methane
	Protonation(Element.Carbon, 1, 1, Hybridization.Sp) to ProtonationTypes("c1", "hc",
		"T", // triple bond
		mapOf(
			Element.Carbon to "c1", // eg ethyne
			Element.Nitrogen to "n1", // TODO hydrogen cyanide?
			Element.Phosphorus to "p2" // TODO methylidynephosphane?
		)
	),
	Protonation(Element.Carbon, 1, 2, Hybridization.Sp2) to ProtonationTypes("c2", "hc",
		"D", // double bond
		mapOf(
			Element.Carbon to "c2", // eg ethene
			Element.Nitrogen to "n2", // TODO ???
			Element.Phosphorus to "p3", // methylenephosphine
			Element.Oxygen to "o" // TODO ???
		)
	),
	Protonation(Element.Carbon, 1, 3, Hybridization.Sp3) to ProtonationTypes("c3", "hc"), // eg ethane
	Protonation(Element.Carbon, 2, 1, Hybridization.Sp2) to ProtonationTypes("c2", "hc"), // eg benzene
	Protonation(Element.Carbon, 2, 2, Hybridization.Sp3) to ProtonationTypes("c3", "hc"), // eg cyclohexane
	Protonation(Element.Carbon, 3, 1, Hybridization.Sp3) to ProtonationTypes("c3", "hc"), // eg isobutane

	Protonation(Element.Nitrogen, 0, 2, Hybridization.Sp2) to ProtonationTypes("n2", "hn"), // eg azanide anion
	Protonation(Element.Nitrogen, 0, 3, Hybridization.Sp3) to ProtonationTypes("n3", "hn"), // eg ammonia
	Protonation(Element.Nitrogen, 0, 4, Hybridization.Sp3) to ProtonationTypes("n4", "hn"), // eg ammonium cation
	Protonation(Element.Nitrogen, 1, 1, Hybridization.Sp) to ProtonationTypes("n1", "hn",
		"T", // triple bond
		mapOf(
			Element.Carbon to "c1", // TODO ???
			Element.Nitrogen to "n1", // eg diazynediium
			Element.Phosphorus to "p2" // TODO ???
		)
	),
	Protonation(Element.Nitrogen, 1, 1, Hybridization.Sp2) to ProtonationTypes("n2", "hn"), // eg diazene
	Protonation(Element.Nitrogen, 1, 2, Hybridization.Sp2) to ProtonationTypes("na", "hn"), // eg formamide
	Protonation(Element.Nitrogen, 1, 2, Hybridization.Sp3) to ProtonationTypes("n3", "hn"), // eg hydrazine
	Protonation(Element.Nitrogen, 1, 3, Hybridization.Sp3) to ProtonationTypes("n4", "hn"), // eg diazanediium
	Protonation(Element.Nitrogen, 2, 1, Hybridization.Sp2) to ProtonationTypes("na", "hn"), // eg N-methylformamide
	Protonation(Element.Nitrogen, 2, 1, Hybridization.Sp3) to ProtonationTypes("n3", "hn"), // eg dimethylamine
	Protonation(Element.Nitrogen, 2, 2, Hybridization.Sp3) to ProtonationTypes("n4", "hn"), // eg dimethylammonium cation
	Protonation(Element.Nitrogen, 3, 1, Hybridization.Sp3) to ProtonationTypes("n4", "hn"), // eg trimethylammonium

	Protonation(Element.Oxygen, 0, 1, Hybridization.Sp2) to ProtonationTypes("o", "hw"), // eg hydroxide anion
	Protonation(Element.Oxygen, 0, 2, Hybridization.Sp3) to ProtonationTypes("ow", "hw"), // eg water (oxidane)
	Protonation(Element.Oxygen, 0, 3, Hybridization.Sp3) to ProtonationTypes("ow", "hw"), // eg hydronium cation
	Protonation(Element.Oxygen, 1, 1, Hybridization.Sp3) to ProtonationTypes("oh", "ho"), // eg hydrogen peroxide
	Protonation(Element.Oxygen, 1, 2, Hybridization.Sp3) to ProtonationTypes("oh", "ho"), // eg methyloxonium cation

	Protonation(Element.Phosphorus, 0, 2, Hybridization.Sp2) to ProtonationTypes("p2", "hp"), // eg phosphanide anion
	Protonation(Element.Phosphorus, 0, 3, Hybridization.Sp3) to ProtonationTypes("p3", "hp"), // eg phosphine
	Protonation(Element.Phosphorus, 0, 4, Hybridization.Sp3) to ProtonationTypes("p5", "hp"), // eg phosphonium cation
	Protonation(Element.Phosphorus, 1, 1, Hybridization.Sp2) to ProtonationTypes("p2", "hp",
		"D", // double bond
		mapOf(
			Element.Carbon to "c2", // eg methylenephosphine
			Element.Nitrogen to "n2", // TODO ???
			Element.Phosphorus to "p3", // TODO diphosphene
			Element.Oxygen to "o" // TODO ???
		)
	),
	/* TODO: can't seem to get PH2+ to be planar here? can't find amber types that work ...
	      let's just not support this protonation for now
	Protonation(Element.Phosphorus, 1, 2, Hybridization.Sp2) to ProtonationTypes("p4", "hp",
		"D", // double bond
		mapOf(
			Element.Carbon to "c2", // eg methylenephosphonium cation
			Element.Nitrogen to "n2", // TODO ???
			Element.Phosphorus to "p3", // TODO ???
			Element.Oxygen to "o" // eg oxophosphonium cation
		)
	),
	*/
	Protonation(Element.Phosphorus, 1, 2, Hybridization.Sp3) to ProtonationTypes("p3", "hp"), // eg diphosphane
	Protonation(Element.Phosphorus, 1, 3, Hybridization.Sp3) to ProtonationTypes("p5", "hp"), // eg methylphosphonium cation
	// TODO: Protonation(Element.Phosphorus, 2, 1, Hybridization.Sp2) ???
	Protonation(Element.Phosphorus, 2, 1, Hybridization.Sp3) to ProtonationTypes("p3", "hp"), // eg dimethylphosphine
	Protonation(Element.Phosphorus, 2, 2, Hybridization.Sp3) to ProtonationTypes("p5", "hp"), // eg dimethylphosphonium
	Protonation(Element.Phosphorus, 3, 1, Hybridization.Sp3) to ProtonationTypes("p5", "hp"), // eg trimethylphosphonium

	Protonation(Element.Sulfur, 0, 1, Hybridization.Sp2) to ProtonationTypes("s", "hs"), // eg bisulfide anion
	Protonation(Element.Sulfur, 0, 2, Hybridization.Sp3) to ProtonationTypes("sh", "hs"), // eg hydrogen sulfide
	Protonation(Element.Sulfur, 0, 3, Hybridization.Sp3) to ProtonationTypes("sh", "hs"), // eg sulfonium cation
	Protonation(Element.Sulfur, 1, 1, Hybridization.Sp3) to ProtonationTypes("sh", "hs"), // eg hydrogen disulfide
	Protonation(Element.Sulfur, 1, 2, Hybridization.Sp3) to ProtonationTypes("sh", "hs") // eg methylsulfonium cation
	// don't think any of the positive oxidation states of sulfur can be protonated, right?
	// wouldn't highly electronegative ligands (eg Oxygen) steal all the protons (ie tautomerization)?
)

enum class Hybridization {
	Sp,
	Sp2,
	Sp3
}
