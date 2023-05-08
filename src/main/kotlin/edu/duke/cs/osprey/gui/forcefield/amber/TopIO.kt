package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.tools.FileTools
import org.joml.Vector3d
import java.util.*
import kotlin.NoSuchElementException
import kotlin.collections.ArrayList
import kotlin.collections.HashMap
import kotlin.math.min


object TopIO {

	/**
	 * Read an AMBER topology file (eg .top)
	 */
	fun read(content: String): AmberTopology {

		val lines = FileTools.parseLines(content)

		// pick out the flags we're interested in
		val flags = HashMap<Flag,FlagInfo>()
		while (lines.hasNext()) {

			val line = lines.next()
			if (line.startsWith("%FLAG ")) {

				// is this a flag we recognize?
				val flag = Flag[line.substring(6).trim()] ?: continue

				// look ahead to the format line
				val fmtline = lines.next().trim()
				val fmt: Any = if (fmtline.startsWith("%FORMAT(")) {
					val fmtspec = fmtline.substring(8, fmtline.length - 1)
					when {
						fmtspec.contains('a') -> StringFormat.of(fmtspec)
						fmtspec.contains('E') -> DoubleFormat.of(fmtspec)
						fmtspec.contains('I') -> IntFormat.of(fmtspec)
						else -> throw ParseException("unrecognized FORMAT: $fmtspec")
					}
				} else {
					throw ParseException("expected FORMAT on line ${lines.lineNumber()}, but got: $fmtline")
				}

				// grab the lines for this flag
				val flagLines = ArrayList<String>()
				while (lines.hasNext() && !lines.peek().startsWith("%")) {
					lines.next()
						.takeUnless { it.isBlank() }
						?.let { flagLines.add(it) }
				}
				flags[flag] = FlagInfo(flagLines, fmt)
			}
		}

		fun get(flag: Flag) = flags[flag] ?: throw NoSuchElementException("missing flag: $flag")

		// parse the flags
		return AmberTopology(
			get(Flag.POINTERS).readInts(),
			get(Flag.ATOM_NAME).readStrings(),
			get(Flag.ATOMIC_NUMBER).readInts(),
			get(Flag.RESIDUE_POINTER).readInts(),
			get(Flag.RESIDUE_LABEL).readStrings(),
			get(Flag.BONDS_INC_HYDROGEN).readInts(),
			get(Flag.BONDS_WITHOUT_HYDROGEN).readInts(),
			get(Flag.DIHEDRALS_INC_HYDROGEN).readInts(),
			get(Flag.DIHEDRALS_WITHOUT_HYDROGEN).readInts(),
			get(Flag.CHARGE).readDoubles(),
			get(Flag.ATOM_TYPE_INDEX).readInts(),
			get(Flag.NONBONDED_PARM_INDEX).readInts(),
			get(Flag.SCEE_SCALE_FACTOR).readDoubles(),
			get(Flag.SCNB_SCALE_FACTOR).readDoubles(),
			get(Flag.LENNARD_JONES_ACOEF).readDoubles(),
			get(Flag.LENNARD_JONES_BCOEF).readDoubles(),
			get(Flag.AMBER_ATOM_TYPE).readStrings()
		)
	}
}

/**
 * A representation of the parts of the AMBER topology file format
 * that we care about in OSPREY.
 *
 * Full docs at: http://ambermd.org/FileFormats.php
 */
class AmberTopology(

	/**
	 * %FORMAT(10i8)  NATOM,    NTYPES, NBONH,  MBONA,  NTHETH, MTHETA,
	 * NPHIH,    MPHIA,  NHPARM, NPARM,  NNB,    NRES,
	 * NBONA,    NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA,
	 * NATYP,    NPHB,   IFPERT, NBPER,  NGPER,  NDPER,
	 * MBPER,    MGPER,  MDPER,  IFBOX,  NMXRS,  IFCAP,
	 * NUMEXTRA, NCOPY
	 */
	val pointers: List<Int>,

	/**
	 * %FORMAT(20a4)  (IGRAPH(i), i=1,NATOM)
	 * IGRAPH : the user-specified atoms names
	 */
	val atomNames: List<String>,

	/**
	 * %FORMAT(10I8)  (ATNUM(i), i=1,NATOM)
	 * ATNUM : the atomic number of each atom.
	 */
	val atomicNumbers: List<Int>,

	/**
	 * %FORMAT(10I8)  (IPRES(i), i=1,NRES)
	 * IPRES  : atoms in each residue are listed for atom "i" in
	 *          IPRES(i) to IPRES(i+1)-1
	 */
	val residuePointers: List<Int>,

	/**
	 * %FORMAT(20A4)  (LBRES(i), i=1,NRES)
	 * LBRES : names of each of the residues
	 */
	val residueLabels: List<String>,

	/**
	 * %FORMAT(10I8)  (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
	 * IBH    : atom involved in bond "i", bond contains hydrogen
	 * JBH    : atom involved in bond "i", bond contains hydrogen
	 * ICBH   : index into parameter arrays RK and REQ
	 */
	val bondsIncHydrogen: List<Int>,

	/**
	 * %FORMAT(10I8)  (IB(i),JB(i),ICB(i), i=1,NBONA)
	 * IB     : atom involved in bond "i", bond does not contain hydrogen
	 * JB     : atom involved in bond "i", bond does not contain hydrogen
	 * ICB    : index into parameter arrays RK and REQ
	 */
	val bondsWithoutHydrogen: List<Int>,

	/**
	 * %FORMAT(10I8)  (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i), i=1,NPHIH)
	 * IPH    : atom involved in dihedral "i", dihedral contains hydrogen
	 * JPH    : atom involved in dihedral "i", dihedral contains hydrogen
	 * KPH    : atom involved in dihedral "i", dihedral contains hydrogen
	 * LPH    : atom involved in dihedral "i", dihedral contains hydrogen
	 * ICPH   : index into parameter arrays PK, PN, PHASE, ONE_SCEE, and ONE_SCNB
	 * for dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i)
	 */
	val dihedralsIncHydrogen: List<Int>,

	/**
	 * %FORMAT(10I8)  (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)
	 * IP     : atom involved in dihedral "i", dihedral does not contain hydrogen
	 * JP     : atom involved in dihedral "i", dihedral does not contain hydrogen
	 * KP     : atom involved in dihedral "i", dihedral does not contain hydrogen
	 * LP     : atom involved in dihedral "i", dihedral does not contain hydrogen
	 * ICP    : index into parameter arrays PK, PN, PHASE, ONE_SCEE, and ONE_SCNB
	 * for dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i).  Note, if the
	 * periodicity is negative, this implies the following entry
	 * in the PK, PN, and PHASE arrays is another term in a
	 * multitermed dihedral.
	 */
	val dihedralsWithoutHydrogen: List<Int>,

	/**
	 * %FORMAT(5E16.8)  (CHARGE(i), i=1,NATOM)
	 * CHARGE : the atom charges.  Amber internally uses units of charge such
	 *          that E = q1*q2/r, where E is in kcal/mol, r is in Angstrom,
	 *          and q1,q2 are the values found in this section of the prmtop file.
	 */
	val charges: List<Double>,

	/**
	 * %FORMAT(1OI8)  (IAC(i), i=1,NATOM)
	 * IAC    : index for the atom types involved in Lennard Jones (6-12)
	 *          interactions.  See ICO below.
	 */
	val atomTypeIndices: List<Int>,

	/**
	 * %FORMAT(10I8)  (ICO(i), i=1,NTYPES*NTYPES)
	 * ICO    : provides the index to the nonbon parameter
	 * arrays CN1, CN2 and ASOL, BSOL.  All possible 6-12
	 * or 10-12 atoms type interactions are represented.
	 * NOTE: A particular atom type can have either a 10-12
	 * or a 6-12 interaction, but not both.  The index is
	 * calculated as follows:
	 * index = ICO(NTYPES*(IAC(i)-1)+IAC(j))
	 * If index is positive, this is an index into the
	 * 6-12 parameter arrays (CN1 and CN2) otherwise it
	 * is an index into the 10-12 parameter arrays (ASOL
	 * and BSOL).
	 */
	val nonbondedParmIndices: List<Int>,

	/**
	 * %FORMAT(5E16.8) (ONE_SCEE(i), i=1,NPTRA)
	 * ONE_SCEE : 1-4 electrostatic scaling constant. It is inverted right after
	 *          it's read in for performance reasons. This allows variable
	 *          1-4 scaling. If not present, it defaults to 1.2 for all dihedrals.
	 *          Therefore, the default ONE_SCEE value in the code is 1.0/1.2
	 */
	val esScaleFactors: List<Double>,

	/**
	 * %FORMAT(5E16.8) (ONE_SCNB(i), i=1,NPTRA)
	 * ONE_SCNB : 1-4 VDW scaling constant. It is inverted right after
	 *          it's read in. This allows variable 1-4 scaling. If not present,
	 *          it defaults to 2.0 for all dihedrals. Therefore, the default
	 *          ONE_SCNB value in the code is 1.0/2.0
	 */
	val vdwScaleFactors: List<Double>,

	/**
	 * %FORMAT(5E16.8)  (CN1(i), i=1,NTYPES*(NTYPES+1)/2)
	 * CN1  : Lennard Jones r**12 terms for all possible atom type interactions,
	 * indexed by ICO and IAC; for atom i and j where i < j, the index
	 * into this array is as follows (assuming the value of ICO(index) is
	 * positive): CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j))).
	 */
	val lennardJonesACoeff: List<Double>,

	/**
	 * %FORMAT(5E16.8)  (CN2(i), i=1,NTYPES*(NTYPES+1)/2)
	 * CN2  : Lennard Jones r**6 terms for all possible
	 * atom type interactions.  Indexed like CN1 above.
	 */
	val lennardJonesBCoeff: List<Double>,

	/**
	 * %FORMAT(20A4)  (ISYMBL(i), i=1,NATOM)
	 * ISYMBL : the AMBER atom types for each atom
	 */
	val atomTypes: List<String>
) {

	init {
		// do some internal consistency checks
		if (residuePointers.size != residueLabels.size) {
			throw IllegalStateException("corrupted topology: inconsistent residue counts: ${residuePointers.size} : ${residueLabels.size}")
		}
	}

	// index the dihedral angle indices so we can find them easily
	class IndexPair(val i1: Int, val i2: Int) {

		override fun hashCode() =
			i1 + i2

		override fun equals(other: Any?) =
			other is IndexPair && (
				(this.i1 == other.i1 && this.i2 == other.i2) // forward order
				|| (this.i1 == other.i2 && this.i2 == other.i1) // reverse order
			)
	}
	val dihedralIndices = HashMap<IndexPair,Int>().apply {
		for (dihedrals in listOf(dihedralsIncHydrogen, dihedralsWithoutHydrogen)) {
			for (i in 0 until dihedrals.size/5) {
				val atomi1 = dihedrals[i*5]
				val atomi2 = dihedrals[i*5 + 3]
				val dihedrali = dihedrals[i*5 + 4]
				this[IndexPair(atomi1, atomi2)] = dihedrali
			}
		}
	}

	fun getDihedralIndex(atomi1: Int, atomi2: Int) =
		dihedralIndices[IndexPair(atomi1, atomi2)]?.let { it - 1 }

	// NOTE: topology file residue "pointers" start counting at 1 instead of 0,
	// so we have to subtract 1 from the "pointer" to get a real index.

	fun firstAtomIndex(resIndex: Int) =
		residuePointers[resIndex] - 1

	fun lastAtomLimit(resIndex: Int) =
		residuePointers.getOrNull(resIndex + 1)?.let { it - 1 } ?: atomNames.size

	fun numAtoms(resIndex: Int) =
		lastAtomLimit(resIndex) - firstAtomIndex(resIndex)

	fun atomIndices(resIndex: Int) =
		firstAtomIndex(resIndex) until lastAtomLimit(resIndex)

	val numAtoms get() = pointers[0]
	val numAtomTypes get() = pointers[1]

	data class MappedAtom(
		val mol: Molecule,
		val res: Polymer.Residue?,
		val atom: Atom?
	)

	fun mapTo(mols: List<Molecule>): Mapped {

		// make a bi-directional mapping between the molecule atoms and the toplogy atoms
		// the topology will have extra atoms, but none should be missing

		/* TODO: this mapping requires that input residues and output residues appear in the same order
		     by default, this seems to be the case, but LEaP does have an option to explicitly
		     preserve the residue order (see reorder_residues). Do we need to turn that on?
		*/

		val molToTop = IdentityHashMap<Atom,Int>()
		val topToMol = HashMap<Int,MappedAtom>()

		var resi = 0

		for (mol in mols) {
			if (mol is Polymer) {

				for (chain in mol.chains) {
					for (res in chain.residues) {

						// get the atom range
						for (i in atomIndices(resi)) {
							val atom = res.atoms.find { it.name == atomNames[i] }
							if (atom != null) {

								// topology atom matched to our molecule
								molToTop[atom] = i
								topToMol[i] = MappedAtom(mol, res, atom)

							} else {

								// topology atom not in our molecule, but still maps to this residue
								topToMol[i] = MappedAtom(mol, null, atom)
							}
						}

						// make sure we matched all the atoms in this residue
						val missingAtoms = res.atoms.filter { it !in molToTop.keys }
						if (missingAtoms.isNotEmpty()) {
							throw IllegalArgumentException("""
								|topology does not describe all atoms in residue ${res.id} ${res.type}:
								|	not found: $missingAtoms
								|	among atoms: ${atomIndices(resi).map { atomNames[it] }}
							""".trimMargin())
						}

						resi += 1
					}
				}

			} else {

				// match this small molecule to the next residue
				for (i in atomIndices(resi)) {
					val atom = mol.atoms.find { it.name == atomNames[i] }
					if (atom != null) {

						// topology atom matched to our molecule
						molToTop[atom] = i
						topToMol[i] = MappedAtom(mol, null, atom)

					} else {

						// topology atom not in our molecule
						topToMol[i] = MappedAtom(mol, null, null)
					}
				}

				// make sure we matched all the atoms in our molecule
				val missingAtoms = mol.atoms.filter { it !in molToTop.keys }
				if (missingAtoms.isNotEmpty()) {
					throw IllegalArgumentException("topology does not describe all atoms in the molecule:\n$missingAtoms")
				}

				resi += 1
			}
		}

		return Mapped(mols, molToTop, topToMol)
	}

	inner class Mapped internal constructor(
		val mols: List<Molecule>,
		private val molToTop: IdentityHashMap<Atom,Int>,
		private val topToMol: HashMap<Int,MappedAtom>
	) {

		/**
		 * Adds all the extra atoms defined by AMBER.
		 *
		 * Returns the number of atoms added.
		 */
		fun addMissingAtoms(
			coords: List<Vector3d>,
			filter: (Atom, Polymer.Residue?) -> Boolean = { _: Atom, _: Polymer.Residue? -> true }
		): Int {

			var numAtomsAdded = 0

			for ((i, mappedAtom) in topToMol) {
				val (mol, res, atom) = mappedAtom

				// skip atoms that already exist in the mol
				if (atom != null) {
					continue
				}

				val newAtom = Atom(
					Element[atomicNumbers[i]],
					atomNames[i],
					Vector3d(coords[i])
				)

				// consult the filter
				if (!filter(newAtom, res)) {
					continue
				}

				// update the molecule (and residue if needed)
				mol.atoms.add(newAtom)
				res?.atoms?.add(newAtom)

				// update the map
				molToTop[newAtom] = i
				topToMol[i] = MappedAtom(mol, res, newAtom)

				numAtomsAdded += 1
			}

			return numAtomsAdded
		}

		/**
		 * Reset all bonds in the molecule to only those defined by the AMBER topology.
		 *
		 * Returns the number of bonds added.
		 */
		fun setBonds(): Int {

			var numBondsAdded = 0

			// combine all the bonds together, we don't need to treat hydrogens specially
			val bonds = bondsIncHydrogen + bondsWithoutHydrogen
			if (bonds.size % 3 != 0) {
				throw IllegalStateException("unexpected number of bond indices: ${bonds.size}, should be a multiple of 3")
			}

			// iterate over the bonds
			val numBonds = bonds.size/3
			for (i in 0 until numBonds) {

				// for some reason, the atom indices are 3x too large. No idea why...
				val i1 = bonds[i*3]/3
				val i2 = bonds[i*3 + 1]/3

				// lookup the atoms for the bond
				// if the atoms haven't been added yet, skip that bond
				val (mol1, _, a1) = topToMol[i1] ?: continue
				val (mol2, _, a2) = topToMol[i2] ?: continue
				a1 ?: continue
				a2 ?: continue

				// make sure the atom names match, just in case
				assert(atomNames[i1] == a1.name)
				assert(atomNames[i2] == a2.name)

				// make sure both atoms are in the same molecule
				if (mol1 !== mol2) {
					throw IllegalStateException("inter-molecular bonds not supported")
				}
				val mol = mol1

				// add the bond
				val wasAdded = mol.bonds.add(a1, a2)
				if (wasAdded) {
					numBondsAdded += 1
				}
			}

			return numBondsAdded
		}

		private fun index(atom: Atom) = molToTop[atom]
			?: throw NoSuchElementException("atom $atom is not in the topology")

		fun atom(i: Int) = topToMol[i]?.atom

		val numAtoms: Int get() = this@AmberTopology.numAtoms

		fun atomType(atom: Atom) = atomTypes[index(atom)]

		fun charge(index: Int) = charges[index]
		fun charge(atom: Atom) = charge(index(atom))

		/**
		 * Returns the electrostatic scaling divisor for the atom pair.
		 * Assumes the atoms are already 1-4 bonded!
		 */
		fun chargeDivisor14(indexa: Int, indexb: Int): Double {
			val i = getDihedralIndex(indexa, indexb) ?: return 1.2
			return esScaleFactors[i]
		}
		fun chargeDivisor14(atoma: Atom, atomb: Atom) = chargeDivisor14(index(atoma), index(atomb))

		fun nonbondIndex(indexa: Int, indexb: Int): Int {
			// ICO(NTYPES*(IAC(i)-1)+IAC(j)) for i < j
			var i = indexa
			var j = indexb
			if (i == j) {
				throw IllegalArgumentException("must have different atoms to get nonbonded index")
			} else if (i > j) {
				val swap = i
				i = j
				j = swap
			}
			val iaci = atomTypeIndices[i]
			val iacj = atomTypeIndices[j]
			// NOTE: Amber uses a lot of 1-offset indexing schemes,
			// so we need to add lots of -1 everywhere to get down to 0-offset indexing.
			// Maybe because amber was based mostly on fortran?
			return nonbondedParmIndices[numAtomTypes*(iaci - 1) + iacj - 1] - 1
		}
		fun nonbondIndex(atoma: Atom, atomb: Atom) = nonbondIndex(index(atoma), index(atomb))

		fun vdw(indexa: Int, indexb: Int): Pair<Double,Double> =
			nonbondIndex(indexa, indexb).let { i ->
				lennardJonesACoeff[i] to lennardJonesBCoeff[i]
			}
		fun vdw(atoma: Atom, atomb: Atom) = vdw(index(atoma), index(atomb))

		/**
		 * Returns the van der Waals scaling divisor for the atom pair.
		 * Assumes the atoms are already 1-4 bonded!
		 */
		fun vdwDivisor14(indexa: Int, indexb: Int): Double {
			val i = getDihedralIndex(indexa, indexb) ?: return 2.0
			return vdwScaleFactors[i]
		}
		fun vdwDivisor14(atoma: Atom, atomb: Atom) = vdwDivisor14(index(atoma), index(atomb))
	}
}

class MappingException(msg: String) : RuntimeException(msg)


private enum class Flag {

	POINTERS,
	ATOM_NAME,
	ATOMIC_NUMBER,
	RESIDUE_POINTER,
	RESIDUE_LABEL,
	BONDS_INC_HYDROGEN,
	BONDS_WITHOUT_HYDROGEN,
	DIHEDRALS_INC_HYDROGEN,
	DIHEDRALS_WITHOUT_HYDROGEN,
	CHARGE,
	ATOM_TYPE_INDEX,
	NONBONDED_PARM_INDEX,
	SCEE_SCALE_FACTOR,
	SCNB_SCALE_FACTOR,
	LENNARD_JONES_ACOEF,
	LENNARD_JONES_BCOEF,
	AMBER_ATOM_TYPE;

	companion object {

		operator fun get(name: String): Flag? =
			values().find { it.name == name }
	}
}

private data class StringFormat(
	val numPerLine: Int,
	val size: Int
) {

	fun parseLine(line: String, into: MutableList<String>) {
		val num = min(line.length/size, numPerLine)
		for (i in 0 until num) {
			into.add(line.substring(i*size, (i+1)*size).trim())
		}
	}

	companion object {
		fun of(spec: String): StringFormat {
			val parts = spec.split('a')
			return StringFormat(
				parts[0].toInt(),
				parts[1].toInt()
			)
		}
	}
}

private data class DoubleFormat(
	val numPerLine: Int,
	val size: Int,
	val precision: Int
) {

	fun parseLine(line: String, into: MutableList<Double>) {
		val num = min(line.length/size, numPerLine)
		for (i in 0 until num) {
			var start = i*size
			val stop = (i+1)*size
			start = line.advancePastWhitespace(start, stop)
			into.add(line.substring(start, stop).trim().toDouble())
		}
	}

	companion object {
		fun of(spec: String): DoubleFormat {
			val parts = spec.split('E', '.')
			return DoubleFormat(
				parts[0].toInt(),
				parts[1].toInt(),
				parts[2].toInt()
			)
		}
	}
}

private data class IntFormat(
	val numPerLine: Int,
	val size: Int
) {
	fun parseLine(line: String, into: MutableList<Int>) {
		val num = min(line.length/size, numPerLine)
		for (i in 0 until num) {
			var start = i*size
			val stop = (i+1)*size
			start = line.advancePastWhitespace(start, stop)
			into.add(line.substring(start, stop).toInt())
		}
	}

	companion object {
		fun of(spec: String): IntFormat {
			val parts = spec.split('I', '.')
			return IntFormat(
				parts[0].toInt(),
				parts[1].toInt()
			)
		}
	}
}

private data class FlagInfo(
	val lines: List<String>,
	val format: Any
) {
	fun readStrings() = ArrayList<String>().apply {
		val format = format as StringFormat
		for (line in lines) {
			format.parseLine(line, this)
		}
	}

	fun readDoubles() = ArrayList<Double>().apply {
		val format = format as DoubleFormat
		for (line in lines) {
			format.parseLine(line, this)
		}
	}

	fun readInts() = ArrayList<Int>().apply {
		val format = format as IntFormat
		for (line in lines) {
			format.parseLine(line, this)
		}
	}
}


class ParseException(msg: String) : RuntimeException(msg)

private fun String.advancePastWhitespace(start: Int, stop: Int): Int {
	var i = start
	while (this[i] == ' ' && start <= stop) {
		i += 1
	}
	return i
}
