package edu.duke.cs.osprey.gui.compiler

import cuchaz.kludge.tools.toRadians
import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.parallelism.Parallelism
import edu.duke.cs.osprey.gui.forcefield.*
import edu.duke.cs.osprey.gui.io.LaunchLimits
import edu.duke.cs.osprey.gui.motions.*
import edu.duke.cs.osprey.gui.prep.Assignments
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.gui.tools.UnsupportedClassException
import edu.duke.cs.osprey.gui.tools.pairs
import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.runBlocking
import org.joml.Vector3d
import java.util.*
import kotlin.collections.ArrayList


/**
 * Collects all the information from a conformation space,
 * combines it with forcefield parameters,
 * and emits a TOML file that describes the parameterized
 * conformation space to Osprey in a totally unambiguous way.
 */
class ConfSpaceCompiler(val confSpace: ConfSpace) {

	val forcefields = ForcefieldSet()
	val netCharges = NetCharges()

	data class Report(
		val warnings: List<CompilerWarning>,
		val error: CompilerError?,
		val compiled: CompiledConfSpace?
	)

	/**
	 * Compiles the conf space and the forcefields.
	 * If something goes wrong, the errors are collected and returned in the compilation report.
	 *
	 * Compilation is actually performed in a separate thread,
	 * but progress can be monitored via the returned progress object.
	 */
	fun compile(): CompilerProgress {

		// get stable orders for all the positions and conformations
		val confSpaceIndex = ConfSpaceIndex(confSpace)

		// make a dummy thread variable for now,
		// so the compiler progress can access it later
		// otherwise we get an intractable dependency cycle
		var thread = null as Thread?

		val numSingles = confSpaceIndex.positions.sumBy { it.fragments.size }
		val numPairs = confSpaceIndex.positions.pairs().sumBy { (posInfo1, posInfo2) ->
			posInfo1.fragments.size*posInfo2.fragments.size
		}

		// init the progress
		val paramsTask = CompilerProgress.Task(
			"Parameterize molecules",
			forcefields.size*(
				confSpaceIndex.mols.size
				+ numSingles
			)
		)
		val fixedAtomsTask = CompilerProgress.Task(
			"Partition fixed atoms",
			forcefields.size*numSingles + 2
		)
		val staticEnergiesTask = CompilerProgress.Task(
			"Calculate energy of static atoms",
			forcefields.size
		)
		val atomPairsTask = CompilerProgress.Task(
			"Calculate forcefield atom pairs",
			1 /* static */ + numSingles + numPairs
		)
		val progress = CompilerProgress(
			paramsTask, fixedAtomsTask, staticEnergiesTask, atomPairsTask,
			threadGetter = { thread ?: throw Error("compilation thread not created yet") }
		)

		// make a thread to do the actual compilation and start it
		thread = Thread {
			compile(confSpaceIndex, progress, paramsTask, fixedAtomsTask, staticEnergiesTask, atomPairsTask)
		}.apply {
			name = "ConfSpaceCompiler"
			isDaemon = false
			start()
		}

		return progress
	}

	private fun compile(
		confSpaceIndex: ConfSpaceIndex,
		progress: CompilerProgress,
		paramsTask: CompilerProgress.Task,
		fixedAtomsTask: CompilerProgress.Task,
		staticEnergiesTask: CompilerProgress.Task,
		atomPairsTask: CompilerProgress.Task
	) {

		val warnings = ArrayList<CompilerWarning>()

		try {

			// limit coroutine max concurrency so we don't run out of memory
			// TODO: make configurable?
			val launchLimits = LaunchLimits(Parallelism.getMaxNumCPUs())

			// compile the forcefield metadata and settings
			val forcefieldInfos = forcefields.map { ff ->
				CompiledConfSpace.ForcefieldInfo(
					ff.forcefield.name,
					ff.forcefield.ospreyImplementation,
					ff.settings().toList()
				)
			}

			// TODO: issues warnings/errors for:
			//   missing forcefield params

			// compute all the forcefield parameters for the conf space
			val params = ConfSpaceParams(confSpaceIndex)

			// first, parameterize the molecule without any conformational changes
			runBlocking {
				for ((moli, mol) in confSpaceIndex.mols.withIndex()) {
					for (ff in forcefields) {
						launchLimits.launch {
							try {
								params[ff, moli] = ff.parameterizeAtoms(
									mol,
									confSpaceIndex.atomIndexWildType(moli),
									netCharges[mol]?.netChargeOrThrow
								)
								paramsTask.increment()
							} catch (t: Throwable) {
								throw CompilerError("Can't parameterize wild-type molecule: $mol", cause = t)
							}
						}
					}
				}
			}

			// then parameterize the molecule for each fragment at each design position
			runBlocking {
				for ((moli, mol) in confSpaceIndex.mols.withIndex()) {
					for (posInfo in confSpaceIndex.positions.filter { it.pos.mol === mol }) {
						posInfo.forEachFrag { assignments, assignmentInfo, confInfo ->
							val atomIndices = confSpaceIndex.matchAtoms(assignments)
							for (ff in forcefields) {
								launchLimits.launch {
									try {
										params[ff, confInfo.fragInfo] = ff.parameterizeAtoms(
											assignmentInfo.molInfo.assignedMol,
											atomIndices[moli],
											netCharges[mol]?.getOrThrow(posInfo.pos, confInfo.fragInfo.frag)
										)
										paramsTask.increment()
									} catch (t: Throwable) {
										throw CompilerError("Can't parameterize molecule: $mol with ${posInfo.pos.name} = ${confInfo.id}", cause = t)
									}
								}
							}
						}
					}
				}
			}

			// analyze all the fixed atoms in the conf space
			val fixedAtoms = FixedAtoms(confSpaceIndex)

			// find the dynamic atoms among the fixed atoms
			// by comparing conf forcefield params against the wild-type forcefield params
			for (ff in forcefields) {
				for (posInfo in confSpaceIndex.positions) {

					val wtParams = params[ff, posInfo.moli]
					val molFixedAtoms = fixedAtoms.fixed(posInfo.moli)

					for (fragInfo in posInfo.fragments) {

						val fragParams = params[ff, fragInfo]

						// find the atoms whose parameters have changed
						val changedAtoms = molFixedAtoms
							.filter { atom ->
								val atomi = confSpaceIndex.atomIndexWildType(posInfo.moli).getOrThrow(atom)
								wtParams[atomi] != fragParams[atomi]
							}

						try {
							fixedAtoms[posInfo].addDynamic(changedAtoms, fragInfo)
						} catch (ex: FixedAtoms.ClaimedAtomException) {

							// convert the exception into a user-friendly compiler error
							val fixedName = ex.atom.fixedName(posInfo.pos.mol)
							val atomi = confSpaceIndex.atomIndexWildType(posInfo.moli)[ex.atom]

							val newDesc = atomi?.let { fragParams[it] } ?: "(unknown params)"

							val oldParams = params[ff, ex.fragInfo]
							val oldDesc = atomi?.let { oldParams[it] } ?: "(unknown params)"

							throw CompilerError("""
								|Forcefield parameterization of fragments at different design positions yielded conflicting parameters for a fixed atom:
								|fixed atom=$fixedName
								|pos1=${ex.fragInfo.posInfo.pos.name}, fragment=${ex.fragInfo.frag.id}, params=$oldDesc
								|pos2=${posInfo.pos.name}, fragment=${fragInfo.frag.id}, params=$newDesc
							""".trimMargin())
						}

						fixedAtomsTask.increment()
					}
				}
			}

			// all the rest of the fixed atoms are static fixed atoms
			fixedAtoms.updateStatic()
			fixedAtomsTask.increment()

			// index the molecules and residues
			val molInfos = confSpaceIndex.mols
				.map { CompiledConfSpace.MolInfo(it.name, it.type) }
			val resIndex = ResidueIndex(confSpaceIndex.mols)

			// compile the static atoms
			val staticAtoms = fixedAtoms.statics
				.map { it.atom.compile(resIndex, it.moli, it.mol) }

			// also map the wild-type atom indices to the static atom indices
			val staticAtomsIndex = HashMap<Int,Int>().apply {
				for (info in fixedAtoms.statics) {
					val wtIndex = confSpaceIndex.atomIndexWildType(info.moli).getOrThrow(info.atom)
					put(wtIndex, info.index)
				}
			}

			fixedAtomsTask.increment()

			// compile the molecule motions
			for ((moli, mol) in confSpaceIndex.mols.withIndex()) {
				for (description in confSpace.molMotions.getOrDefault(mol, mutableListOf())) {
					molInfos[moli].motions.add(description.compile(mol, fixedAtoms))
				}
			}

			// collect the atoms affected by molecule motions
			val staticAffectedAtomsByMol = confSpaceIndex.mols
				.map { mol ->
					val used = Atom.identitySet()
					confSpace.molMotions[mol]
						?.flatMap { it.getAffectedAtoms() }
						?.filter { atom ->

							// intersect with the static atoms
							// (some molecule motions can move dynamic atoms too)
							if (!fixedAtoms.isStatic(atom)) {
								return@filter false
							}

							// de-duplicate
							if (atom in used) {
								return@filter false
							}
							used.add(atom)

							true
						}
						?: emptyList()
				}

			// calculate the internal energies for the static atoms that aren't affected by molecule motions
			val staticUnaffectedAtomsByMol = fixedAtoms.staticAtomsByMol
				.mapIndexed { moli, atoms ->
					atoms.filter { it !in staticAffectedAtomsByMol[moli] }
				}
			val staticEnergies = runBlocking {
				forcefields
					.map { ff ->
						async {

							// start with the internal energies for affected atoms
							var staticEnergy = staticAffectedAtomsByMol
								.mapIndexed { moli, atoms ->
									val molParams = params[ff, moli]
									val atomIndices = confSpaceIndex.atomIndexWildType(moli)
									atoms.mapNotNull { atom ->
										val atomi = atomIndices.getOrThrow(atom)
										molParams[atomi]?.internalEnergy()
									}.sum()
								}.sum()

							// add the pairwise energies for unaffected atoms
							staticEnergy += ForcefieldCalculator.calc(
								ff.parameterizeAtomPairs(staticUnaffectedAtomsByMol.mapIndexed { moli, _ ->
									ForcefieldParams.MolInfo(moli, confSpaceIndex.mols[moli], params[ff, moli], confSpaceIndex.atomIndexWildType(moli))
								}),
								staticUnaffectedAtomsByMol.mapIndexed { moli, atoms ->
									ForcefieldCalculator.MolInfo(moli, confSpaceIndex.mols[moli], atoms, confSpaceIndex.atomIndexWildType(moli), params[ff, moli])
								},
								ff.unconnectedDistance
							)

							return@async staticEnergy
								.also { staticEnergiesTask.increment() }
						}
					}
					.awaitAll()
			}

			// compile the design positions
			val posInfos = ArrayList<CompiledConfSpace.PosInfo>()
			runBlocking {
				for (posInfo in confSpaceIndex.positions) {

					// compile the fragments
					val fragInfos = ArrayList<CompiledConfSpace.FragInfo>()
					posInfo.forEachFrag { assignments, assignmentInfo, confInfo ->
						launchLimits.launch {
							fragInfos.add(confInfo.fragInfo.compile(fixedAtoms, assignmentInfo))
						}
					}

					// compile the conformations
					val confInfos = ArrayList<CompiledConfSpace.ConfInfo>()
					posInfo.forEachConf { assignments, assignmentInfo, confInfo ->
						launchLimits.launch {
							confInfos.add(confInfo.compile(resIndex, assignments, assignmentInfo, confInfo, confSpaceIndex, fixedAtoms, params))
						}
					}

					posInfos.add(CompiledConfSpace.PosInfo(
						posInfo.pos.name,
						posInfo.pos.type,
						fragInfos,
						confInfos
					))
				}
			}

			// compile all the atom pairs for the forcefields
			val atomPairs = forcefields
				.map { AtomPairs(confSpaceIndex) }
			runBlocking {

				// compile the static-self pairs
				launchLimits.launch {
					atomPairs.compileStatic(confSpaceIndex, staticAtomsIndex, fixedAtoms, staticAffectedAtomsByMol, staticUnaffectedAtomsByMol, params)
					atomPairsTask.increment()
				}

				for (posInfo1 in confSpaceIndex.positions) {

					posInfo1.forEachFrag { assignments, assignmentInfo1, confInfo1 ->

						// compile the pos-self and pos-static atom pairs
						launchLimits.launch {
							atomPairs.compilePos(assignments, confInfo1, assignmentInfo1, confSpaceIndex, staticAtomsIndex, fixedAtoms, params)
							atomPairsTask.increment()
						}
					}

					for (posInfo2 in confSpaceIndex.positions.subList(0, posInfo1.index)) {
						confSpace.forEachFragIn(posInfo1, posInfo2) { assignments, assignmentInfo1, confInfo1, assignmentInfo2, confInfo2 ->

							// compile the pos-pos atom pairs
							launchLimits.launch {
								atomPairs.compilePosPos(assignments, confInfo1, assignmentInfo1, confInfo2, assignmentInfo2, confSpaceIndex, fixedAtoms, params)
								atomPairsTask.increment()
							}
						}
					}
				}
			}

			// if we made it this far, return a successful compiler result
			progress.report = Report(
				warnings,
				null,
				CompiledConfSpace(
					confSpace.name,
					forcefieldInfos,
					molInfos,
					resIndex.resInfos,
					staticAtoms,
					staticEnergies,
					posInfos,
					atomPairs
				)
			)

		} catch (t: Throwable) {

			// dump the error to the console, just in case any developers are watching
			t.printStackTrace(System.err)

			// wrap the exception in a compiler error if needed
			val error = t as? CompilerError
				?: CompilerError("Error", null, t)

			// return a report with only warnings and errors, but no compiled conf space
			progress.report = Report(warnings, error, null)
		}
	}

	private class ResidueIndex(mols: List<Molecule>) {

		// index all the residues
		val resInfos: List<CompiledConfSpace.ResInfo>
		val resIndices: Map<String,Int>
		init {
			val resInfos = ArrayList<CompiledConfSpace.ResInfo>()
			resIndices = HashMap()
			for (mol in mols.filterIsInstance<Polymer>()) {
				for (chain in mol.chains) {
					for ((resi, res) in chain.residues.withIndex()) {
						resIndices[res.id] = resInfos.size
						resInfos.add(CompiledConfSpace.ResInfo(chain.id, res.id, res.type, resi))
					}
				}
			}
			this.resInfos = resInfos
		}

		fun indexOfRes(res: Polymer.Residue?): Int =
			if (res != null) {
				resIndices[res.id] ?: throw IllegalArgumentException("residue has no index: $res")
			} else {
				-1
			}
	}

	/**
	 * Compiles the motions for this molecule
	 */
	private fun MolMotion.Description.compile(
		mol: Molecule,
		fixedAtoms: FixedAtoms
	): CompiledConfSpace.MotionInfo = when (this) {

		is DihedralAngle.MolDescription -> CompiledConfSpace.MotionInfo.DihedralAngle(
			minDegrees,
			maxDegrees,
			abcd = listOf(a, b, c, d).map { fixedAtoms.getStatic(it).index },
			rotated = rotatedAtoms.map { fixedAtoms.getStatic(it).index }
		)

		is TranslationRotation.MolDescription -> CompiledConfSpace.MotionInfo.TranslationRotation(
			maxTranslationDist,
			maxRotationDegrees.toRadians(),
			centroid = Vector3d().apply {
				for (atom in mol.atoms) {
					add(atom.pos)
				}
				div(mol.atoms.size.toDouble())
			}
		)

		else -> throw UnsupportedClassException("don't know how to compile molecule motion", this)
	}


	/**
	 * Compiles the fragment
	 */
	private fun ConfSpaceIndex.FragInfo.compile(fixedAtoms: FixedAtoms, assignmentInfo: Assignments.AssignmentInfo): CompiledConfSpace.FragInfo {

		// collect the atoms for this fragment in order, including dynamic fixed atoms
		val confAtoms = orderAtoms(fixedAtoms, assignmentInfo)

		return CompiledConfSpace.FragInfo(
			name = frag.id,
			atomNames = confAtoms.map { it.name }
		)
	}

	/**
	 * Compiles the conformation
	 */
	private fun ConfSpaceIndex.ConfInfo.compile(
		resIndex: ResidueIndex,
		assignments: Assignments,
		assignmentInfo: Assignments.AssignmentInfo,
		confInfo: ConfSpaceIndex.ConfInfo,
		confSpaceIndex: ConfSpaceIndex,
		fixedAtoms: FixedAtoms,
		params: ConfSpaceParams
	): CompiledConfSpace.ConfInfo {

		// collect the atoms for this fragment in order, including dynamic fixed atoms
		val confAtoms = fragInfo.orderAtoms(fixedAtoms, assignmentInfo)

		// compile the motions for this conformation
		val confAtomsIndex = AtomIndex(confAtoms)
		val motions = confConfSpace.motions
			.map { it.compile(confAtomsIndex, fixedAtoms, assignmentInfo) }

		// compute the internal energies of the conformation atoms
		val atomIndices = confSpaceIndex.matchAtoms(assignments)
		val internalEnergies = forcefields.map { ff ->
			val fragParams = params[ff, fragInfo]
			confAtoms.mapNotNull { atom ->
				val atomi = atomIndices[confInfo.posInfo.moli].getOrThrow(atom)
				fragParams[atomi]?.internalEnergy()
			}
			.sum()
		}

		return CompiledConfSpace.ConfInfo(
			id = id,
			type = fragInfo.frag.type,
			atoms = confAtoms.map { it.compile(resIndex, posInfo.moli, assignmentInfo.molInfo.assignedMol) },
			fragIndex = fragInfo.index,
			motions = motions,
			internalEnergies = internalEnergies
		)
	}

	/**
	 * Compiles the continuous motion
	 */
	private fun ConfMotion.Description.compile(
		confAtoms: AtomIndex,
		fixedAtoms: FixedAtoms,
		assignmentInfo: Assignments.AssignmentInfo
	): CompiledConfSpace.MotionInfo {

		// get the index of the static atom (but make it negative, in a reversible way)
		fun FixedAtoms.getStaticNegated(atom: Atom, assignmentInfo: Assignments.AssignmentInfo): Int =
			getStatic(assignmentInfo.molInfo.getConfSpaceAtomOrThrow(atom)).let { -it.index - 1 }

		when (this) {

			is DihedralAngle.ConfDescription -> {

				val a = assignmentInfo.confSwitcher.atomResolverOrThrow.resolveOrThrow(motion.a)
				val b = assignmentInfo.confSwitcher.atomResolverOrThrow.resolveOrThrow(motion.b)
				val c = assignmentInfo.confSwitcher.atomResolverOrThrow.resolveOrThrow(motion.c)
				val d = assignmentInfo.confSwitcher.atomResolverOrThrow.resolveOrThrow(motion.d)

				return CompiledConfSpace.MotionInfo.DihedralAngle(
					minDegrees,
					maxDegrees,
					abcd = listOf(a, b, c, d)
						.map { confAtoms[it] ?: fixedAtoms.getStaticNegated(it, assignmentInfo) },
					rotated = DihedralAngle.findRotatedAtoms(assignmentInfo.molInfo.assignedMol, b, c)
						.map { confAtoms.getOrThrow(it) }
				)
			}

			else -> throw UnsupportedClassException("don't know how to compile conformation motion", this)
		}
	}

	private fun Atom.compile(resIndex: ResidueIndex, moli: Int, mol: Molecule): CompiledConfSpace.AtomInfo =
		CompiledConfSpace.AtomInfo(
			name,
			Vector3d(pos.checkForErrors(name)),
			moli,
			resIndex.indexOfRes((mol as? Polymer)?.findResidue(this))
		)

	/**
	 * Compile atom pairs for pos-self and pos-static interactions.
	 */
	private suspend fun List<AtomPairs>.compilePos(
		assignments: Assignments,
		confInfo: ConfSpaceIndex.ConfInfo,
		assignmentInfo: Assignments.AssignmentInfo,
		confSpaceIndex: ConfSpaceIndex,
		staticAtomsIndex: Map<Int,Int>,
		fixedAtoms: FixedAtoms,
		params: ConfSpaceParams
	) {
		val atomIndices = confSpaceIndex.matchAtoms(assignments)

		// parameterize
		val atomPairsParams = forcefields.map { ff ->
			// OPTIMIZATION: this *should be* the slow part
			ff.parameterizeAtomPairs(confSpaceIndex.mols.withIndex().map { (moli, mol) ->
				val molInfo = assignments.molInfoByConfSpaceMol(mol)
				val atomsParams = if (molInfo === assignmentInfo.molInfo) {
					// use fragment-specific params for this molecule, since this design position has been assigned
					params[ff, confInfo.fragInfo]
				} else {
					// otherwise, use the wild-type params for this molecule
					params[ff, moli]
				}
				ForcefieldParams.MolInfo(moli, molInfo.assignedMol, atomsParams, atomIndices[moli])
			})
		}

		// position side
		val posAtoms = confInfo.fragInfo.orderAtoms(fixedAtoms, assignmentInfo)
		val posAtomIndex = AtomIndex(posAtoms)
		val posMolInfos = listOf(
			AtomPairer.MolInfo(confInfo.posInfo.moli, assignmentInfo.molInfo.assignedMol, posAtoms, atomIndices[confInfo.posInfo.moli])
		)

		// static side
		val staticAtoms = fixedAtoms.staticAtomsByMol(assignments)
		val staticMolInfos = confSpaceIndex.mols.withIndex().map { (moli, mol) ->
			AtomPairer.MolInfo(moli, assignments.molInfoByConfSpaceMol(mol).assignedMol, staticAtoms[moli], atomIndices[moli])
		}

		// pos-self interactions
		for (molPair in AtomPairer.molPairs(posMolInfos, posMolInfos)) {
			molPair.forEach(forcefields.unconnectedDistance()) { info1, atomi1, info2, atomi2, distance ->

				// map the global atom indices to fragment-local atom indices
				val atomiPos1 = posAtomIndex.getOrThrow(info1.atomIndex.getOrThrow(atomi1))
				val atomiPos2 = posAtomIndex.getOrThrow(info2.atomIndex.getOrThrow(atomi2))

				for (ffi in forcefields.indices) {
					val atomPairParams = atomPairsParams[ffi][info1.moli, atomi1, info2.moli, atomi2, distance] ?: continue
					this[ffi].pos.add(
						confInfo.posInfo.index,
						confInfo.fragInfo.index,
						atomiPos1,
						atomiPos2,
						atomPairParams.list
					)
				}
			}
		}

		// pos-static interactions
		for (molPair in AtomPairer.molPairs(posMolInfos, staticMolInfos)) {
			molPair.forEach(forcefields.unconnectedDistance()) { infoPos, atomiPos, infoStatic, atomiWildType, distance ->

				// map the global atom indices to fragment-local atom indices
				val atomiFrag = posAtomIndex.getOrThrow(infoPos.atomIndex.getOrThrow(atomiPos))

				// map the static atom's index from the wild-type set to the static set
				val atomiStatic  = staticAtomsIndex.getStaticIOrThrow(atomiWildType)

				for (ffi in forcefields.indices) {
					val atomPairParams = atomPairsParams[ffi][infoPos.moli, atomiPos, infoStatic.moli, atomiWildType, distance] ?: continue
					this[ffi].posStatic.add(
						confInfo.posInfo.index,
						confInfo.fragInfo.index,
						atomiFrag,
						atomiStatic,
						atomPairParams.list
					)
				}
			}
		}
	}

	/**
	 * Compile atom pairs for pos-pos interactions.
	 */
	private suspend fun List<AtomPairs>.compilePosPos(
		assignments: Assignments,
		confInfo1: ConfSpaceIndex.ConfInfo,
		assignmentInfo1: Assignments.AssignmentInfo,
		confInfo2: ConfSpaceIndex.ConfInfo,
		assignmentInfo2: Assignments.AssignmentInfo,
		confSpaceIndex: ConfSpaceIndex,
		fixedAtoms: FixedAtoms,
		params: ConfSpaceParams
	) {
		val atomIndices = confSpaceIndex.matchAtoms(assignments)

		// position 1
		val confAtoms1 = confInfo1.fragInfo.orderAtoms(fixedAtoms, assignmentInfo1)
		val confAtomIndex1 = AtomIndex(confAtoms1)
		val molInfos1 = listOf(
			AtomPairer.MolInfo(confInfo1.posInfo.moli, assignmentInfo1.molInfo.assignedMol, confAtoms1, atomIndices[confInfo1.posInfo.moli])
		)

		// position 2
		val confAtoms2 = confInfo2.fragInfo.orderAtoms(fixedAtoms, assignmentInfo2)
		val confAtomIndex2 = AtomIndex(confAtoms2)
		val molInfos2 = listOf(
			AtomPairer.MolInfo(confInfo2.posInfo.moli, assignmentInfo2.molInfo.assignedMol, confAtoms2, atomIndices[confInfo2.posInfo.moli])
		)

		// parameterize
		val atomPairsParams = if (confInfo1.posInfo.moli == confInfo2.posInfo.moli) {

			// design positions are in the same molecule, need to get combined atoms params
			val moli = confInfo1.posInfo.moli
			val molInfo = assignmentInfo1.molInfo
			val atomIndex = atomIndices[confInfo1.posInfo.moli]

			fun List<FixedAtoms.DynamicInfo>.dynamicIndices(): Set<Int> =
				map { atomIndex.getOrThrow(molInfo.getAssignedAtomOrThrow(it.atom)) }
				.toSet()

			fun Set<Atom>.wtConfIndices(): Set<Int> =
				map { confSpaceIndex.atomIndexWildType(moli).getOrThrow(it) }
				.toSet()

			forcefields.map { ff ->
				// combine the atoms params for both confs together
				// use the dynamic atoms to break ties in atom params
				// and ignore the wild-type conf atoms from the other positions
				val atomsParams = ff.combineAtomsParams(
					ForcefieldParams.AtomsInfo(params[ff, confInfo1.fragInfo], fixedAtoms[confInfo1.posInfo].dynamics.dynamicIndices(), confInfo2.posInfo.pos.sourceAtoms.wtConfIndices()),
					ForcefieldParams.AtomsInfo(params[ff, confInfo2.fragInfo], fixedAtoms[confInfo2.posInfo].dynamics.dynamicIndices(), confInfo1.posInfo.pos.sourceAtoms.wtConfIndices())
				)

				// OPTIMIZATION: this *should be* the slow part
				ff.parameterizeAtomPairs(listOf(
					ForcefieldParams.MolInfo(moli, molInfo.assignedMol, atomsParams, atomIndex)
				))
			}

		} else {

			// design positions are in different molecules, can just use cached atoms params
			forcefields.map { ff ->
				// OPTIMIZATION: this *should be* the slow part
				ff.parameterizeAtomPairs(listOf(
					ForcefieldParams.MolInfo(confInfo1.posInfo.moli, assignmentInfo1.molInfo.assignedMol, params[ff, confInfo1.fragInfo], atomIndices[confInfo1.posInfo.moli]),
					ForcefieldParams.MolInfo(confInfo2.posInfo.moli, assignmentInfo2.molInfo.assignedMol, params[ff, confInfo2.fragInfo], atomIndices[confInfo2.posInfo.moli])
				))
			}
		}

		// for each atom pair ...
		for (molPair in AtomPairer.molPairs(molInfos1, molInfos2)) {
			molPair.forEach(forcefields.unconnectedDistance()) { info1, atomi1, info2, atomi2, distance ->

				// map the global atom indices to fragment-local atom indices
				val atomiConf1 = confAtomIndex1.getOrThrow(info1.atomIndex.getOrThrow(atomi1))
				val atomiConf2 = confAtomIndex2.getOrThrow(info2.atomIndex.getOrThrow(atomi2))

				for (ffi in forcefields.indices) {
					val atomPairParams = atomPairsParams[ffi][info1.moli, atomi1, info2.moli, atomi2, distance] ?: continue
					this[ffi].posPos.add(
						confInfo1.posInfo.index,
						confInfo1.fragInfo.index,
						confInfo2.posInfo.index,
						confInfo2.fragInfo.index,
						atomiConf1,
						atomiConf2,
						atomPairParams.list
					)
				}
			}
		}
	}

	/**
	 * Compile atom pairs for static self-interactions
	 */
	private suspend fun List<AtomPairs>.compileStatic(
		confSpaceIndex: ConfSpaceIndex,
		staticAtomsIndex: Map<Int,Int>,
		fixedAtoms: FixedAtoms,
		staticAffectedAtomsByMol: List<List<Atom>>,
		staticUnaffectedAtomsByMol: List<List<Atom>>,
		params: ConfSpaceParams
	) {
		val wildTypeAtomsIndices = confSpaceIndex.mols.indices.map { confSpaceIndex.atomIndexWildType(it) }
		val affectedMolInfos = fixedAtoms.mols.mapIndexed { moli, mol ->
			AtomPairer.MolInfo(moli, mol, staticAffectedAtomsByMol[moli], wildTypeAtomsIndices[moli])
		}
		val unaffectedMolInfos = fixedAtoms.mols.mapIndexed { moli, mol ->
			AtomPairer.MolInfo(moli, mol, staticUnaffectedAtomsByMol[moli], wildTypeAtomsIndices[moli])
		}

		// do affected-self interactions, and affected-unaffected interactions
		val groups = listOf(
			affectedMolInfos to affectedMolInfos,
			affectedMolInfos to unaffectedMolInfos
		)

		// OPTIMIZATION: this *should be* the slow part
		val atomPairsParams = forcefields.map { ff ->
			ff.parameterizeAtomPairs(fixedAtoms.mols.mapIndexed { moli, mol ->
				ForcefieldParams.MolInfo(moli, mol, params[ff, moli], wildTypeAtomsIndices[moli])
			})
		}

		// for each atom pair ...
		for ((molInfos1, molInfos2) in groups) {

			for (molPair in AtomPairer.molPairs(molInfos1, molInfos2)) {
				molPair.forEach(forcefields.unconnectedDistance()) { info1, atomi1, info2, atomi2, distance ->

					// map the static atom's index from the wild-type set to the static set
					val atomiStatic1  = staticAtomsIndex.getStaticIOrThrow(atomi1)
					val atomiStatic2  = staticAtomsIndex.getStaticIOrThrow(atomi2)

					for (ffi in forcefields.indices) {
						val atomPairParams = atomPairsParams[ffi][info1.moli, atomi1, info2.moli, atomi2, distance] ?: continue
						this[ffi].static.add(
							atomiStatic1,
							atomiStatic2,
							atomPairParams.list
						)
					}
				}
			}
		}
	}
}

private fun Map<Int,Int>.getStaticIOrThrow(atomiWildType: Int): Int =
	this[atomiWildType] ?: throw NoSuchElementException("no static index for wild-type atom $atomiWildType")
