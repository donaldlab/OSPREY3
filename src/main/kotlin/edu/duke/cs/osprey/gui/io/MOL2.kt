package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.molscope.molecule.*
import org.joml.Vector3d
import java.text.DecimalFormat
import java.text.FieldPosition
import java.util.*
import kotlin.NoSuchElementException
import kotlin.collections.ArrayList
import kotlin.collections.HashMap


// this file handles IO with the Mol2 format:
// http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf


class Mol2Metadata {
	val atomTypes: MutableMap<Atom,String> = IdentityHashMap()
	val atomCharges: MutableMap<Atom,String> = IdentityHashMap()
	val bondTypes: MutableMap<AtomPair,String> = HashMap()
	val dictionaryTypes: MutableMap<Polymer.Residue,String> = IdentityHashMap()
	var smallMoleculeDictionaryType: String = defaultSmallMoleculeDictionaryType

	companion object {

		const val defaultCharge = "0.0"
		const val defaultBondType = "1" // 1 is single bond
		const val defaultDictionaryType = "1" // 1 is protein
		const val defaultSmallMoleculeDictionaryType = "0" // ??? no idea, hopefully it's nothing bad
	}
}


/**
 * Save the molecule to the Mol2 format.
 */
fun Molecule.toMol2(metadata: Mol2Metadata? = null): String {

	val mol = this

	val buf = StringBuffer()
	fun write(str: String, vararg args: Any) = buf.append(String.format(str, *args))

	// get the chains (either from the polymer, or if a small molecule, make a dummy chain)
	val smallMoleculeRes: Polymer.Residue?
	val chains: List<Polymer.Chain>
	if (mol is Polymer) {
		smallMoleculeRes = null
		chains = mol.chains
	} else {
		smallMoleculeRes = Polymer.Residue("1", mol.type ?: "MOL", mol.atoms)
		chains = listOf(Polymer.Chain("A").apply {
			residues.add(smallMoleculeRes)
		})
	}

	// assign an id to each residue
	val resIds = IdentityHashMap<Polymer.Residue,Int>().apply {
		var i = 1
		for (chain in chains) {
			for (res in chain.residues) {
				put(res, i)
				i += 1
			}
		}
	}

	// get the residue for each atom
	val residuesByAtom = IdentityHashMap<Atom,Polymer.Residue>().apply {
		for (chain in chains) {
			for (res in chain.residues) {
				for (atom in res.atoms) {
					put(atom, res)
				}
			}
		}
	}
	fun Atom.getResidue() =
		residuesByAtom[this] ?: throw NoSuchElementException("No residue found for atom $this in $mol")

	// get all the atom indices (1-based)
	val indicesByAtom = IdentityHashMap<Atom,Int>()
	mol.atoms.forEachIndexed { i, atom ->
		indicesByAtom[atom] = i + 1
	}

	// flatten the bonds to a list
	val bonds = mol.bonds.toSet().toList()

	// write the molecule section
	write("@<TRIPOS>MOLECULE\n")
	write(mol.name.replace("\n", " ") + "\n")
	write("  ${mol.atoms.size} ${bonds.size} ${chains.sumBy { it.residues.size }}\n")
	when (this) {
		is Polymer -> write("BIOPOLYMER\n")
		else -> write("SMALL\n")
	}
	write("NO_CHARGES\n")
	write("\n")
	write("\n")

	val types = metadata?.atomTypes
	fun Atom.getType() =
		if (types != null) {
			types[this] ?: throw NoSuchElementException("No atom type found for atom $this @ ${residuesByAtom[this]} in $mol")
		} else {
			// no types explicitly defined for this molecule, so probably that means no one cares what we write here
			// but we have to write something, so just use the atom element
			element.symbol
		}

	val charges = metadata?.atomCharges
	fun Atom.getCharge() =
		if (charges != null) {
			charges[this] ?: throw NoSuchElementException("No charge found for atom $this @ ${residuesByAtom[this]} in $mol")
		} else {
			// no charges explicitly defined for this molecule, so probably that means no one cares what we write here
			// but we have to write something, so just use an arbitrary value
			Mol2Metadata.defaultCharge
		}

	// write out atom coords with up to 6 digits of precision
	val coordsFormat = DecimalFormat("#.000000")

	// when this reflection voodoo actually works, use it to get to the fast versions of the Double->String functions
	// tragically, an oversight in the DecimalFormat API prevents us from using the fast functions via normal means  T_T
	// yes, profiling showed the original slow functions (String.format) was actually a noticeable performance bottleneck
	val dontCareFieldPosition = try {
		Class.forName("java.text.DontCareFieldPosition").run {
			getDeclaredField("INSTANCE").run {
				isAccessible = true
				get(null) as FieldPosition
			}
		}
	} catch (t: Throwable) {
		null
	}
	fun StringBuffer.appendCoord(value: Double) =
		if (dontCareFieldPosition != null) {
			coordsFormat.format(value, this, dontCareFieldPosition)
		} else {
			append(coordsFormat.format(value))
		}

	// write the atom section
	write("@<TRIPOS>ATOM\n")
	for (atom in atoms) {
		val res = atom.getResidue()
		buf.append(' ')
		buf.append(' ')
		buf.append(indicesByAtom[atom])
		buf.append(' ')
		buf.append(atom.name)
		buf.append(' ')
		buf.appendCoord(atom.pos.x)
		buf.append(' ')
		buf.appendCoord(atom.pos.y)
		buf.append(' ')
		buf.appendCoord(atom.pos.z)
		buf.append(' ')
		buf.append(atom.getType())
		buf.append(' ')
		buf.append(resIds[res])
		buf.append(' ')
		buf.append(res.type)
		buf.append(' ')
		buf.append(atom.getCharge())
		buf.append('\n')
	}

	// write the bond section
	write("@<TRIPOS>BOND\n")
	bonds.forEachIndexed { i, bond ->
		buf.append(' ')
		buf.append(' ')
		buf.append(i + 1)
		buf.append(' ')
		buf.append(indicesByAtom.getValue(bond.a))
		buf.append(' ')
		buf.append(indicesByAtom.getValue(bond.b))
		buf.append(' ')
		buf.append(metadata?.bondTypes?.getValue(bond) ?: Mol2Metadata.defaultBondType)
		buf.append('\n')
	}

	// write the substructure section
	write("@<TRIPOS>SUBSTRUCTURE\n")
	for (chain in chains) {
		for (res in chain.residues) {
			buf.append(' ')
			buf.append(' ')
			buf.append(resIds[res])
			buf.append(' ')
			buf.append(res.id)
			buf.append(' ')
			buf.append(indicesByAtom[res.atoms.first()])
			buf.append(" RESIDUE ")
			buf.append(if (res == smallMoleculeRes) {
				metadata?.smallMoleculeDictionaryType ?: Mol2Metadata.defaultSmallMoleculeDictionaryType
			} else {
				metadata?.dictionaryTypes?.getValue(res) ?: Mol2Metadata.defaultDictionaryType
			})
			buf.append(' ')
			buf.append(chain.id)
			buf.append(' ')
			buf.append(res.type)
			buf.append('\n')
		}
	}

	return buf.toString()
}

/**
 * Read a molecule in MOL2 format.
 */
fun Molecule.Companion.fromMol2(mol2: String): Molecule = fromMol2WithMetadata(mol2).first

fun Molecule.Companion.fromMol2WithMetadata(mol2: String): Pair<Molecule,Mol2Metadata> {

	// parse a few sections from the mol2 file
	val lines = mol2.lines()
	val sections = HashMap<Section,List<String>>()
	var i = 0
	fun peek() = lines.getOrNull(i) ?: ""
	fun next() = peek().also { i += 1 }
	while (i < lines.size) {

		val line = next()
		if (line.startsWith("@")) {

			// is this a section we care about?
			val section = Section[
				line.trim()
					.substring(1)
					.toUpperCase()
					.replace("<", "")
					.replace(">", "_")
			] ?: continue

			// grab the lines for this section
			val sectionLines = ArrayList<String>()
			while (i < lines.size && !peek().startsWith("@")) {
				next()
					.takeUnless { it.isBlank() }
					?.let { sectionLines.add(it) }
			}
			sections[section] = sectionLines
		}
	}

	val metadata = Mol2Metadata()

	// build a molecule from the sections
	val sectionMol = sections[Section.TRIPOS_MOLECULE]
		?: throw NoSuchElementException("missing MOLECULE section")
	val molName = sectionMol[0].trim()

	val sectionSub = sections[Section.TRIPOS_SUBSTRUCTURE]
		?: throw NoSuchElementException("missing SUBSTRUCTURE section")

	// if there's more than one substructure, assume the molecule is a polymer
	val mol = if (sectionSub.size > 1) {
		Polymer(molName)
	} else {

		// try to get the mol type from the (single) substructure
		val parts = sectionSub[0].tokenize()

		Molecule(molName, parts[6])
	}

	// read the atoms
	val sectionAtoms = sections[Section.TRIPOS_ATOM]
		?: throw NoSuchElementException("missing ATOM section")
	val atomsById = HashMap<String,Atom>()
	val atomsBySub = HashMap<String,ArrayList<Atom>>()
	for (line in sectionAtoms) {
		val parts = line.tokenize()

		val id = parts[0]
		val name = parts[1]
		val pos = Vector3d(
			parts[2].toDoubleOrNull() ?: throw Mol2ParseException(line, "\"${parts[2]}\" doesn't seem to be a number"),
			parts[3].toDoubleOrNull() ?: throw Mol2ParseException(line, "\"${parts[3]}\" doesn't seem to be a number"),
			parts[4].toDoubleOrNull() ?: throw Mol2ParseException(line, "\"${parts[4]}\" doesn't seem to be a number")
		)
		val type = parts[5]
		val subId = parts.getOrNull(6)
		val charge = parts.getOrNull(8)

		val element = Element.findByPrefixMatch(type)
			?: Element.findByPrefixMatch(name)
			?: throw NoSuchElementException("can't determine element for atom $name $type")
			// haha, this exception seems aptly named =P

		val atom = Atom(element, name, pos)
		atomsById[id] = atom
		mol.atoms.add(atom)

		// track the substructure, if needed
		if (subId != null) {
			atomsBySub.computeIfAbsent(subId) { ArrayList() }.add(atom)
		}

		// update metadata
		metadata.atomTypes[atom] = type
		if (charge != null) {
			metadata.atomCharges[atom] = charge
		}
	}

	// read the bonds
	val sectionBonds = sections[Section.TRIPOS_BOND]
		?: throw NoSuchElementException("missing BOND section")

	for (line in sectionBonds) {
		val parts = line.tokenize()

		val i1 = parts[1]
		val i2 = parts[2]
		val type = parts[3]

		val a1 = atomsById[i1] ?: throw NoSuchElementException("no atom with id $i1")
		val a2 = atomsById[i2] ?: throw NoSuchElementException("no atom with id $i2")

		val bond = AtomPair(a1, a2)
		mol.bonds.add(bond)

		// update the metadata
		metadata.bondTypes[bond] = type
	}

	// read the polmyer, if any
	if (mol is Polymer) {

		for (line in sectionSub) {
			val parts = line.tokenize()

			val id = parts[0]
			val name = parts[1]
			val dictionaryType = parts.getOrNull(4)
			val chainId = parts.getOrNull(5) ?: "A"
			val subType = parts.getOrNull(6) ?: "MOL"

			val chain = mol.chains
				.find { it.id == chainId }
				?: Polymer.Chain(chainId).apply {
					mol.chains.add(this)
				}

			val atoms = atomsBySub[id] ?: ArrayList()

			val res = Polymer.Residue(name, subType, atoms)
			chain.residues.add(res)

			// update the metadata
			if (dictionaryType != null) {
				metadata.dictionaryTypes[res] = dictionaryType
			}
		}
	}

	return mol to metadata
}


private fun String.tokenize() = split(" ").filter { it.isNotBlank() }


private enum class Section {

	TRIPOS_MOLECULE,
	TRIPOS_ATOM,
	TRIPOS_BOND,
	TRIPOS_SUBSTRUCTURE;

	companion object {

		operator fun get(name: String): Section? =
			values().find { it.name == name }
	}
}


class Mol2ParseException(line: String, msg: String) : RuntimeException("$msg\nin line:\n$line")
