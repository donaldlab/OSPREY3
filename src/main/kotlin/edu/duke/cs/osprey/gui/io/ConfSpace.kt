package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.forcefield.amber.findTypeOrThrow
import edu.duke.cs.osprey.gui.motions.DihedralAngle
import edu.duke.cs.osprey.gui.motions.TranslationRotation
import edu.duke.cs.osprey.gui.prep.Anchor
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.gui.prep.DesignPosition
import edu.duke.cs.osprey.gui.tools.UnsupportedClassException
import org.tomlj.Toml
import org.tomlj.TomlPosition
import java.util.*
import kotlin.NoSuchElementException
import kotlin.collections.ArrayList
import kotlin.collections.HashMap


private const val wtConflibId = "__wildtype__"

fun ConfSpace.toToml(): String {

	val buf = StringBuilder()
	fun write(str: String, vararg args: Any) = buf.append(String.format(str, *args))

	// get the molecules in a stable order
	val mols = mols.map { (_, mol) -> mol }

	// write out the molecules, and keep the atom indices
	val (molsToml, indicesByAtom) = mols.toOMOLMapped(flattenAtomIndices = true)
	write(molsToml)

	fun Atom.index() =
		indicesByAtom[this]
			?: throw NoSuchElementException("no index for atom $this")

	val fragConflibIds = IdentityHashMap<ConfLib.Fragment,String>()

	// make a conflib for the wild-type conformations
	val wtConflib = ConfLib(
		id = wtConflibId,
		name = "Wild-Type",
		fragments = wildTypeFragments().associateBy { it.id }
	)
	for (frag in wtConflib.fragments.values) {
		fragConflibIds[frag] = wtConflib.id
	}

	// write down all the conf libs
	for (conflib in conflibs + listOf(wtConflib)) {

		for (frag in conflib.fragments.values) {
			fragConflibIds[frag] = conflib.id
		}

		write("\n")
		write(conflib.toToml(table = "conflibs.${conflib.id}"))
	}

	fun ConfLib.Fragment.conflibId() =
		fragConflibIds[this]
			?: throw Error("no conformation library id for fragment: $id\n" +
				"maybe the conformation library was not added to the conformation space")

	// write the conf space
	write("\n")
	write("[confspace]\n")
	write("name = %s\n", name.quote())

	val posIndices = HashMap<DesignPosition,Int>()
	for ((moli, mol) in mols.withIndex()) {

		// write the design positions, if any
		designPositionsByMol[mol]?.forEach pos@{ pos ->

			// get the position conf space, or skip this pos
			val posConfSpace = positionConfSpaces[pos] ?: return@pos

			val posi = posIndices.size
			posIndices[pos] = posi

			write("\n")
			write("[confspace.positions.$posi]\n")
			write("mol = %s\n", moli)
			write("name = %s\n", pos.name.quote())
			write("type = %s\n", pos.type.quote())

			// write the atoms
			write("atoms = [ %s ]\n",
				pos.sourceAtoms
					.map { it.index() }
					.sorted()
					.joinToString(", ")
			)

			// write the anchors
			for ((i, group) in pos.anchorGroups.withIndex()) {
				write("\n")
				write("[confspace.positions.$posi.confspace.anchorGroups.$i]\n")
				write("anchors = [\n")
				for (anchor in group) {
					val type = when (anchor) {
						is Anchor.Single -> "single"
						is Anchor.Double -> "double"
					}
					write("\t{ type = %s, atoms = [ %s ] },\n",
						type.quote(),
						anchor.anchorAtoms
							.map { it.index() }
							.joinToString(", ")
					)
				}
				write("]\n")
			}

			// write the position conf space
			write("\n")
			write("[confspace.positions.$posi.confspace]\n")
			posConfSpace.wildTypeFragment?.let {
				write("wtfrag = %s\n", it.id.quote())
			}
			write("mutations = [ %s ]\n",
				posConfSpace.mutations.joinToString(", ") { it.quote() }
			)

			// write the conf conf spaces (in a sorted order)
			for ((confi, space) in posConfSpace.confs.withIndex()) {

				write("[confspace.positions.$posi.confspace.confs.$confi]\n")
				write("conflib = %s\n", space.frag.conflibId().quote())
				write("frag = %s\n", space.frag.id.quote())
				write("conf = %s\n", space.conf.id.quote())

				// write the motions, if any
				if (space.motions.isNotEmpty()) {
					write("motions = [\n")
					for (motion in space.motions) {
						when (motion) {

							is DihedralAngle.ConfDescription -> {
								write("\t{ type = %s, index = %d, minDegrees = %12.6f, maxDegrees = %12.6f },\n",
									"dihedralAngle".quote(),
									space.frag.motions
										.indexOf(motion.motion)
										.takeIf { it >= 0 }
										?: throw NoSuchElementException("dihedral angle motion is not in fragment"),
									motion.minDegrees,
									motion.maxDegrees
								)
							}

							else -> throw UnsupportedClassException("don't know how to save conformation motion", motion)
						}
					}
					write("]\n")
				}
			}
		}

		// write the molecule motions, if any
		molMotions[mol]?.forEachIndexed { motioni, motion ->
			write("\n")
			write("[confspace.molMotions.$moli.$motioni]\n")
			when (motion) {

				is DihedralAngle.MolDescription -> {
					write("type = %s\n", "dihedralAngle".quote())
					write("atoms = [ %s ]\n",
						listOf(motion.a, motion.b, motion.c, motion.d)
							.joinToString(", ") { it.index().toString() }
					)
					write("degrees = [ %12.6f, %12.6f ]\n", motion.minDegrees, motion.maxDegrees)
				}

				is TranslationRotation.MolDescription -> {
					write("type = %s\n", "translationRotation".quote())
					write("maxTranslationDist = %12.6f\n", motion.maxTranslationDist)
					write("maxRotationDegrees = %12.6f\n", motion.maxRotationDegrees)
				}

				else -> throw UnsupportedClassException("don't know how to save molecule motion", motion)
			}
		}
	}

	return buf.toString()
}

private fun String.quote() = "'$this'"

fun ConfSpace.Companion.fromToml(toml: String): ConfSpace {

	// parse the TOML
	val doc = Toml.parse(toml)
	if (doc.hasErrors()) {
		throw TomlParseException("TOML parsing failure:\n${doc.errors().joinToString("\n")}")
	}

	// read the molecules
	val molsAndAtoms = Molecule.fromOMOLWithAtoms(doc)
	fun getMol(i: Int, pos: TomlPosition?) =
		molsAndAtoms.getOrNull(i)
			?.let { (mol, _) -> mol }
			?: throw TomlParseException("no molecule found with index $i", pos)

	// calculate the types for the molecules
	val typesAndMols =
		molsAndAtoms
			.map { (mol, _) -> mol.findTypeOrThrow() to mol }

	// build the atom lookup
	val atomIndices = HashMap<Int,Atom>().apply {
		for ((_, atomIndices) in molsAndAtoms) {
			for ((i, atom) in atomIndices) {
				put(i, atom)
			}
		}
	}
	fun getAtom(i: Int, pos: TomlPosition?) =
		atomIndices[i]
			?: throw TomlParseException("no atom with index $i", pos)

	// read the conf space
	val confSpaceTable = doc.getTableOrThrow("confspace")

	val confSpace = ConfSpace(typesAndMols).apply {

		name = confSpaceTable.getStringOrThrow("name")

		// read the conflibs, if any
		var wtConflib: ConfLib? = null
		val conflibsTable = doc.getTable("conflibs")
		if (conflibsTable != null) {
			for (conflibId in conflibsTable.keySet()) {
				val conflibTable = conflibsTable.getTableOrThrow(conflibId)

				val conflib = ConfLib.from(conflibTable)

				// add all but the wildtype library to the conf space
				if (conflibId == wtConflibId) {
					wtConflib = conflib
				} else {
					conflibs.add(conflib)
				}
			}
		}

		// read the positions, if any
		val positionsTable = confSpaceTable.getTable("positions")
		if (positionsTable != null) {
			for (poskey in positionsTable.keySet()) {
				val posTable = positionsTable.getTableOrThrow(poskey)
				val posPos = positionsTable.inputPositionOf(poskey)

				val moli = posTable.getIntOrThrow("mol")
				val mol = getMol(moli, posPos)

				designPositionsByMol.getOrPut(mol) { ArrayList() }.add(DesignPosition(
					name = posTable.getStringOrThrow("name", posPos),
					type = posTable.getStringOrThrow("type", posPos),
					mol = mol
				).apply pos@{

					// read atoms
					val atomsArray = posTable.getArrayOrThrow("atoms", posPos)
					for (i in 0 until atomsArray.size()) {
						sourceAtoms.add(getAtom(atomsArray.getInt(i), posPos))
					}

					// read the pos conf space
					val posConfSpaceTable = posTable.getTableOrThrow("confspace", posPos)
					val posConfSpacePos = posTable.inputPositionOf("confspace")

					positionConfSpaces.getOrMake(this@pos).apply {

						// get the wild type fragment, if any
						wildTypeFragment = posConfSpaceTable.getString("wtfrag")
							?.let { id -> wtConflib?.fragments?.get(id) }

						// read the mutations
						val mutationsArray = posConfSpaceTable.getArrayOrThrow("mutations", posConfSpacePos)
						for (i in 0 until mutationsArray.size()) {
							mutations.add(mutationsArray.getString(i))
						}

						// read the confs
						val confsTable = posConfSpaceTable.getTable("confs")
						if (confsTable != null) {
							for (confkey in confsTable.keySet()) {
								val confTable = confsTable.getTableOrThrow(confkey)
								val confPos = confsTable.inputPositionOf(confkey)

								// get the conformation library
								val conflibId = confTable.getStringOrThrow("conflib", confPos)
								val conflib =
									if (conflibId == wtConflibId) {
										wtConflib ?: throw TomlParseException("Wild-type fragment was specified, but no wild-type conformation library was included")
									} else {
										conflibs.getOrThrow(conflibId)
									}

								// get the fragment
								val fragId = confTable.getStringOrThrow("frag", confPos)
								val frag = conflib.fragments[fragId]
									?: throw TomlParseException("no fragment with id $fragId in conformation library $conflibId")

								// get the conf
								val confId = confTable.getStringOrThrow("conf")
								val conf = frag.confs[confId]
									?: throw TomlParseException("no conf with id $confId in fragment ${frag.id}", confPos)

								// make the conf conf space
								confs.add(frag, conf).apply {

									// read the motions, if any
									val motionsArray = confTable.getArray("motions")
									if (motionsArray != null) {
										for (motioni in 0 until motionsArray.size()) {
											val motionTable = motionsArray.getTable(motioni)
											val motionPos = motionsArray.inputPositionOf(motioni)

											when (motionTable.getStringOrThrow("type", motionPos)) {

												"dihedralAngle" -> {
													val index = motionTable.getIntOrThrow("index", motionPos)
													motions.add(DihedralAngle.ConfDescription(
														this@pos,
														(frag.motions.getOrNull(index) as? ConfLib.ContinuousMotion.DihedralAngle)
															?: throw NoSuchElementException("no dihedral motion at fragment $frag at index $index"),
														motionTable.getDoubleOrThrow("minDegrees"),
														motionTable.getDoubleOrThrow("maxDegrees")
													))
												}
											}
										}
									}
								}
							}
						}
					}

					// read anchor groups
					val anchorGroupsTable = posConfSpaceTable.getTableOrThrow("anchorGroups", posConfSpacePos)
					for (anchorGroupKey in anchorGroupsTable.keySet()) {
						val anchorGroupTable = anchorGroupsTable.getTableOrThrow(anchorGroupKey)
						val anchorGroupPos = anchorGroupsTable.inputPositionOf(anchorGroupKey)
						val anchorArray = anchorGroupTable.getArrayOrThrow("anchors", anchorGroupPos)

						anchorGroups.add(ArrayList<Anchor>().apply {

							for (i in 0 until anchorArray.size()) {
								val anchorTable = anchorArray.getTable(i)
								val anchorPos = anchorArray.inputPositionOf(i)

								val anchorAtomsArray = anchorTable.getArrayOrThrow("atoms", anchorPos)

								when (anchorTable.getStringOrThrow("type", anchorPos)) {
									"single" ->
										add(anchorSingle(
											a = getAtom(anchorAtomsArray.getInt(0), anchorPos),
											b = getAtom(anchorAtomsArray.getInt(1), anchorPos),
											c = getAtom(anchorAtomsArray.getInt(2), anchorPos)
										))
									"double" ->
										add(anchorDouble(
											a = getAtom(anchorAtomsArray.getInt(0), anchorPos),
											b = getAtom(anchorAtomsArray.getInt(1), anchorPos),
											c = getAtom(anchorAtomsArray.getInt(2), anchorPos),
											d = getAtom(anchorAtomsArray.getInt(3), anchorPos)
										))
								}
							}
						})
					}
				})
			}
		}

		// read the molecule motions, if any
		val molMotionsTable = confSpaceTable.getTable("molMotions")
		if (molMotionsTable != null) {

			for (molkey in molMotionsTable.keySet()) {
				val moli = molkey.toInt()
				val (_, mol) = mols[moli]

				val motionsTable = molMotionsTable.getTable(molkey)
				if (motionsTable != null) {

					for (motionkey in motionsTable.keySet()) {
						val motionTable = motionsTable.getTableOrThrow(motionkey)
						val motionPos = motionsTable.inputPositionOf(motionkey)

						val motions = molMotions.getOrPut(mol) { ArrayList() }

						when (motionTable.getString("type")) {

							"dihedralAngle" -> {

								val atomsArray = motionTable.getArrayOrThrow("atoms", motionPos)
								val degreesArray = motionTable.getArrayOrThrow("degrees", motionPos)

								motions.add(DihedralAngle.MolDescription(
									mol,
									getAtom(atomsArray.getInt(0), motionPos),
									getAtom(atomsArray.getInt(1), motionPos),
									getAtom(atomsArray.getInt(2), motionPos),
									getAtom(atomsArray.getInt(3), motionPos),
									degreesArray.getDouble(0),
									degreesArray.getDouble(1)
								))
							}

							"translationRotation" -> {
								motions.add(TranslationRotation.MolDescription(
									mol,
									motionTable.getDoubleOrThrow("maxTranslationDist"),
									motionTable.getDoubleOrThrow("maxRotationDegrees")
								))
							}
						}
					}
				}
			}
		}
	}

	return confSpace
}
