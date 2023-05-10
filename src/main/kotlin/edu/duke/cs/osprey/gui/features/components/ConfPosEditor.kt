package edu.duke.cs.osprey.gui.features.components

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.columns
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.gui.infoTip
import edu.duke.cs.osprey.molscope.render.MoleculeRenderEffects
import edu.duke.cs.osprey.molscope.view.MoleculeRenderStack
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.features.slide.ConformationEditor.MolInfo
import edu.duke.cs.osprey.gui.features.slide.ConformationEditor.PosInfo
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.motions.DihedralAngle
import edu.duke.cs.osprey.gui.prep.Assignments
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.gui.prep.PosAssignment
import edu.duke.cs.osprey.gui.prep.assign
import kotlin.random.asKotlinRandom


class ConfPosEditor(val confSpace: ConfSpace, val molInfo: MolInfo, val posInfo: PosInfo, val onClose: () -> Unit) {

	private val winState = WindowState()
		.apply { pOpen.value = true }
	private val posEditor = DesignPositionEditor(confSpace, posInfo.pos)
	private var resetTabSelection = true
	private val discreteTabState = Commands.TabState()
	private val continuousTabState = Commands.TabState()

	private val pDihedralRadiusDegrees = Ref.of(
		getter = { posEditor.dihedralSettings.radiusDegrees.toFloat() },
		setter = { posEditor.dihedralSettings.radiusDegrees = it.toDouble() }
	)
	private val pDihedralIncludeHydroxyls = Ref.of(posEditor.dihedralSettings::includeHydroxyls)
	private val pDihedralIncludeNonHydroxylHGroups = Ref.of(posEditor.dihedralSettings::includeNonHydroxylHGroups)

    private inner class MutInfo(val type: String) {

        inner class FragInfo(val conflib: ConfLib?, val frag: ConfLib.Fragment) {

            val id =
                if (conflib != null) {
                    frag.uniqueId(conflib)
                } else {
                    "wildtype/${frag.id}"
                }

            inner class ConfInfo(val conf: ConfLib.Conf) {

                val fragInfo get() = this@FragInfo

                val id =
                    if (conflib != null) {
                        conf.uniqueId(conflib, frag)
                    } else {
                        "wildtype/${frag.id}/${conf.id}"
                    }
                val label = "${conf.name}##$id"

                val isSelected = posInfo.posConfSpace.confs.contains(frag, conf)

                fun add() {
                    posInfo.posConfSpace.confs.add(frag, conf).run {

                        // add default dihedral angles from the conformation library
                        DihedralAngle.ConfDescription.makeFromLibrary(
                            posInfo.pos,
                            frag,
                            conf,
                            posEditor.dihedralSettings
                        )
                            .forEach { motions.add(it) }
                    }
                }

                fun remove() = posInfo.posConfSpace.confs.remove(frag, conf)
                fun get() = posInfo.posConfSpace.confs.get(frag, conf)
                fun included() = posInfo.posConfSpace.confs.contains(frag, conf)
            }

            // collect the conformations
            val confInfos = frag.confs
                .values
                .map { ConfInfo(it) }

            val numSelected = confInfos.count { it.isSelected }
            val numPossible = frag.confs.size

            val label = "${frag.name}, $numSelected/$numPossible confs###${posInfo.pos.name}-$id"
        }

        // collect the fragments from the relevant conflib
        val fragInfos = confSpace.getConflibsByMol(molInfo.mol)
            .flatMap { conflib -> conflib.fragments.values.map { conflib to it } }
            .toMutableList()
            .let { frags ->
                // prepend the wild-type fragment, if any
                val wtFrag = posInfo.posConfSpace.wildTypeFragment
                if (wtFrag != null) {
                    listOf(null to wtFrag) + frags
                } else {
                    frags
                }
            }
            .filter { (_, frag) -> frag.type == type && posInfo.pos.isFragmentCompatible(frag) }
            .map { (conflib, frag) -> FragInfo(conflib, frag) }

        val numSelected = fragInfos.sumOf { it.numSelected }
        val numConfs = fragInfos.sumOf { it.numPossible }

        val label = "$type, $numSelected/$numConfs confs###${posInfo.pos.name}-$type"
    }
	private val mutInfos = ArrayList<MutInfo>()
	private var selectedFragInfo: MutInfo.FragInfo? = null
	private var selectedConfInfo: MutInfo.FragInfo.ConfInfo? = null

	private val motionInfos = ArrayList<MotionInfo>()
	private val pJiggle = Ref.of(false)

	private val rand = java.util.Random().asKotlinRandom()

	private var stackedMol: MoleculeRenderStack.StackedMol? = null
	private var renderEffects = null as MoleculeRenderEffects.Writer?

	fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		val view = slide.views
			.filterIsInstance<MoleculeRenderView>()
			.find { it.molStack.originalMol == molInfo.mol }
			?: throw Error("can't init design position, molecule has no render view")

		winState.render(
			onOpen = {

				// init the pos editor
				posEditor.init(view)
			},
			whenOpen = {

				// draw the window
				window("Edit ${posInfo.pos.name} conformations###${slide.name}_confPosEditor_${posInfo.label}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					// if we need to reset the tab selection, make the flags for the first tab
					fun makeFlags() =
						if (resetTabSelection) {
							resetTabSelection = false
							IntFlags.of(Commands.TabItemFlags.SetSelected)
						} else {
							IntFlags.of(Commands.TabItemFlags.None)
						}

					tabBar("tabs") {

						// only show the position editing gui for flexible positions
						if (posInfo.isFlexible) {
							when (molInfo.molType) {
								MoleculeType.Protein -> tabItem("Protein", flags = makeFlags()) {
									posEditor.guiProtein(imgui, slidewin, view)
								}
								// TODO: others?
								else -> Unit
							}
							tabItem("Atoms", flags = makeFlags()) {
								posEditor.guiAtoms(imgui, slidewin, view)
							}
						}

						tabItem(discreteTabState, "Discrete", flags = makeFlags(),
							onActivated = {
								activateDiscreteTab()
							},
							whenActive = {
								renderDiscreteTab(imgui, view)
							},
							onDeactivated = {
								deactivateDiscreteTab(view)
							}
						)

						tabItem(continuousTabState, "Continuous", flags = makeFlags(),
							onActivated = {
								activateContinuousTab(view)
							},
							whenActive = {
								renderContinuousTab(imgui, view)
							},
							onDeactivated = {
								deactivateContinuousTab(view)
							}
						)
					}
				}
			},
			onClose = {

				// deactivate any tabs if they're active
				if (discreteTabState.wasActive) {
					discreteTabState.wasActive = false
					deactivateDiscreteTab(view)
				}
				if (continuousTabState.wasActive) {
					continuousTabState.wasActive = false
					deactivateContinuousTab(view)
				}

				// remove any molecule and effects render overrides
				stackedMol?.pop()
				stackedMol = null
				renderEffects?.close()
				renderEffects = null

				// cleanup the pos editor
				posEditor.closed()
				onClose()
			}
		)
	}

	private fun resetInfos() {
		mutInfos.clear()
		for (mut in posInfo.posConfSpace.mutations) {
			mutInfos.add(MutInfo(mut))
		}
	}

	private fun activateDiscreteTab() {

		resetInfos()

		// find the wildtype conf info, if possible, and select it by default
		val wtConf = posInfo.posConfSpace.wildTypeFragment
			?.confs
			?.values
			?.firstOrNull()
		selectedConfInfo = mutInfos
			.flatMap { mutInfo -> mutInfo.fragInfos.flatMap { it.confInfos } }
			.find { it.conf === wtConf }
	}

	private fun renderDiscreteTab(imgui: Commands, view: MoleculeRenderView) = imgui.run {

		var changed = false

		text("Conformations:")
		sameLine()
		infoTip("""
				|Select the conformations you want to include in your design by clicking the checkboxes.
				|You can temporarily preview a conformation by selecting the radion button next to it.
				|All temporary confomrations will be reverted when you're finished with the conformation editor.
			""".trimMargin())
		child("confs", 300f, 400f, true) {
			for (mutInfo in mutInfos) {
				treeNode(mutInfo.label, flags = IntFlags.of(Commands.TreeNodeFlags.DefaultOpen)) {
					for (fragInfo in mutInfo.fragInfos) {
						treeNode(fragInfo.label) {
							for (confInfo in fragInfo.confInfos) {

								if (radioButton("##radio-${confInfo.id}", selectedConfInfo === confInfo)) {
									selectedConfInfo = confInfo
									setConf(view, fragInfo.frag, confInfo.conf)
								}
								sameLine()
								checkbox(confInfo.label, confInfo.isSelected)?.let { isChecked ->

									// mark the conf as included or not in the design position conf space
									if (isChecked) {
										confInfo.add()
									} else {
										confInfo.remove()
									}

									changed = true
								}
							}
						}
					}
				}
			}
			if (mutInfos.isEmpty()) {
				text("(No compatible conformations)")
			}
		}

		if (button("Select all")) {
			for (mutInfo in mutInfos) {
				for (fragInfo in mutInfo.fragInfos) {
					for (confInfo in fragInfo.confInfos) {
						if (!confInfo.included()) {
							confInfo.add()
						}
					}
				}
			}
			changed = true
		}
		sameLine()
		if (button("Deselect all")) {
			for (mutInfo in mutInfos) {
				for (fragInfo in mutInfo.fragInfos) {
					for (confInfo in fragInfo.confInfos) {
						confInfo.remove()
					}
				}
			}
			changed = true
		}

		// apply changes if needed
		if (changed) {
			resetInfos()
		}
	}

	private fun deactivateDiscreteTab(@Suppress("UNUSED_PARAMETER") view: MoleculeRenderView) {

		// remove any molecule and effects render overrides
		stackedMol?.pop()
		stackedMol = null
		renderEffects?.close()
		renderEffects = null

		// reset the pos editor to re-draw the render effects
		posEditor.resetInfos()
	}

	private fun setConf(view: MoleculeRenderView, frag: ConfLib.Fragment, conf: ConfLib.Conf) {

		// make the assignment
		val posAssignment = PosAssignment(posInfo.pos, frag, conf)
		val assignmentInfo = confSpace.assign(posAssignment).assignmentInfos.getValue(posAssignment)
		val mol = assignmentInfo.molInfo.assignedMol

		// override the molecule render view with the new molecule
		stackedMol?.pop()
		stackedMol = view.molStack.push(mol)
		if (renderEffects == null) {
			renderEffects = view.renderEffects.writer()
		}

		// add render effects to match the position editor's style
		assignmentInfo.confSwitcher.anchorMatch?.posAnchors?.forEach { anchor ->
			for (atom in anchor.anchorAtoms) {
				renderEffects?.set(atom, DesignPositionEditor.anchorEffect)
			}
		}
		for (atom in assignmentInfo.confSwitcher.currentAtoms) {
			renderEffects?.set(atom, DesignPositionEditor.selectedEffect)
		}
	}

	private fun selectConf(view: MoleculeRenderView, fragInfo: MutInfo.FragInfo, confInfo: MutInfo.FragInfo.ConfInfo) {

		selectedFragInfo = fragInfo
		selectedConfInfo = confInfo

		setConf(view, fragInfo.frag, confInfo.conf)

		resetMotionInfos(view)
	}

	private fun activateContinuousTab(view: MoleculeRenderView) {

		resetInfos()

		// find the wildtype fragment, if possible, and select it by default
		val wtFragInfo = posInfo.posConfSpace.wildTypeFragment?.let { frag ->
			for (mutInfo in mutInfos) {
				for (fragInfo in mutInfo.fragInfos) {
					if (fragInfo.frag == frag) {
						return@let fragInfo
					}
				}
			}
			return@let null
		}
		val wtConfInfo = wtFragInfo
			?.confInfos
			?.firstOrNull()
		if (wtFragInfo != null && wtConfInfo != null) {
			selectConf(view, wtFragInfo, wtConfInfo)
		} else {
			selectedFragInfo = null
			selectedConfInfo = null
			resetMotionInfos(view)
		}
	}

	private fun resetMotionInfos(view: MoleculeRenderView) {

		motionInfos.clear()

		// get the selected fragment,conformation
		val fragInfo = selectedFragInfo ?: return
		val confInfo = selectedConfInfo ?: return
		val space = confInfo.get() ?: return

		// make an assignment with this conf, so we can modify coordinates without changing the original molecules
		val assignment = PosAssignment(posInfo.pos, fragInfo.frag, confInfo.conf)
		val assignmentInfo = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

		// override the molecule render view with the new molecule
		stackedMol?.pop()
		stackedMol = view.molStack.push(assignmentInfo.molInfo.assignedMol)
		if (renderEffects == null) {
			renderEffects = view.renderEffects.writer()
		}

		// add render effects to match the position editor's style
		assignmentInfo.confSwitcher.anchorMatch?.posAnchors?.forEach { anchor ->
			for (atom in anchor.anchorAtoms) {
				renderEffects?.set(atom, DesignPositionEditor.anchorEffect)
			}
		}
		for (atom in assignmentInfo.confSwitcher.currentAtoms) {
			renderEffects?.set(atom, DesignPositionEditor.selectedEffect)
		}

		for (motion in fragInfo.frag.motions) {
			when (motion) {

				is ConfLib.ContinuousMotion.DihedralAngle -> {
					motionInfos.add(MotionInfo.DihedralInfo(
						view,
						motion,
						posInfo,
						assignmentInfo,
						space,
						posEditor.dihedralSettings.radiusDegrees
					))
				}
			}
		}
	}

	private fun renderContinuousTab(imgui: Commands, view: MoleculeRenderView) = imgui.run {

		// HACKHACK: ImGUI doesn't want to resize the window automatically for this tab for some reason
		// so put a prop of a certain size to push out the width of the window
		child("foo", 600f, 1f) {}

		columns(2) {
			column(300f) { // frags/confs column

				text("Choose a fragment:")
				child("testFrags", 280f, 270f, true) {

					var hasSelections = false

					for (mutInfo in mutInfos) {
						for (fragInfo in mutInfo.fragInfos) {

							// only show fragments with selected confs
							val firstConfInfo = fragInfo.confInfos
								.firstOrNull { it.isSelected }
								?: continue

							if (radioButton("${fragInfo.frag.name}##fragRadio-${fragInfo.id}", selectedFragInfo === fragInfo)) {
								selectConf(view, fragInfo, firstConfInfo)
							}

							hasSelections = true
						}
					}

					if (!hasSelections) {
						text("No fragments to show here.")
						text("Try including some conformations")
						text("at this design position.")
					}
				}

				text("Choose a conformation:")
				child("testConfs", 280f, 270f, true) {
					val fragInfo = selectedFragInfo
					if (fragInfo != null) {
						for (confInfo in fragInfo.confInfos) {

							// only show selected confs
							if (!confInfo.included()) {
								continue
							}

							if (radioButton("${confInfo.conf.name}##confRadio-${confInfo.id}", selectedConfInfo === confInfo)) {
								selectConf(view, fragInfo, confInfo)
							}
						}
					}
				}
			}
			column(300f) col@{ // motions column

				val confInfo = selectedConfInfo ?: return@col

				text("Dihedral Angle Settings:")
				indent(20f)

				text("Angle radius (in degrees)")
				sameLine()
				infoTip("""
						|This value specifies the half-width (ie radius) of the interval of allowed angles.
					""".trimMargin())

				sliderFloat("##dihedralRadius", pDihedralRadiusDegrees, 0f, 180f, "%.1f")
				sameLine()
				infoTip("Ctrl-click to type a precise value")

				spacing()

				checkbox("Include Hydroxyls", pDihedralIncludeHydroxyls)

				spacing()

				checkbox("Include Non-Hydroxyl H-groups", pDihedralIncludeNonHydroxylHGroups)
				sameLine()
				infoTip("e.g. Methyls, Methylenes")

				fun ConfSpace.ConfConfSpace.resetDihedralAngles() {

					// recreate all the dihedral angles with the new settings
					motions.removeIf { it is DihedralAngle.ConfDescription }
					motions.addAll(DihedralAngle.ConfDescription.makeFromLibrary(
						posInfo.pos,
						frag,
						conf,
						posEditor.dihedralSettings
					))
				}

				if (button("Reset ${confInfo.fragInfo.frag.name}/${confInfo.conf.name}")) {
					confInfo.get()?.resetDihedralAngles()
					resetMotionInfos(view)
				}
				sameLine()
				infoTip("Resets dihedral angles for just this selected conformation.")

				if (button("Reset ${confInfo.fragInfo.frag.name}/*")) {
					confInfo.fragInfo.confInfos
						.mapNotNull { it.get() }
						.forEach { it.resetDihedralAngles() }
					resetMotionInfos(view)
				}
				sameLine()
				infoTip("Resets dihedral angles for all conformations in this selected fragment.")

				if (button("Reset */*")) {
					mutInfos
						.flatMap { it.fragInfos }
						.flatMap { it.confInfos }
						.mapNotNull { it.get() }
						.forEach { it.resetDihedralAngles() }
					resetMotionInfos(view)
				}
				sameLine()
				infoTip("Resets dihedral angles for all conformations in this design position.")


				unindent(20f)

				spacing()

				text("Degrees of Freedom for: ${confInfo.fragInfo.frag.name}/${confInfo.conf.name}")
				indent(20f)

				if (motionInfos.isNotEmpty()) {

					for ((i, motionInfo) in motionInfos.withIndex()) {

						// breathe a little
						spacing()
						spacing()
						separator()
						spacing()
						spacing()

						withId(i) {
							text(motionInfo.label)
							motionInfo.gui(imgui)
						}
					}

					spacing()
					spacing()
					separator()
					spacing()
					spacing()

					// this is my new favorite button! =D
					checkbox("Jiggle", pJiggle)
					sameLine()
					infoTip("""
							|Randomly sample all the degrees of freedom over time,
							|to get a sense of the the range of conformations accessible
							|to the design position.
						""".trimMargin())

					if (pJiggle.value) {
						motionInfos.forEach { it.viewer?.jiggle(rand) }
					}

				} else {

					text("(No degrees of freedom found)")
				}

				unindent(20f)
			}
		}
	}

	private fun deactivateContinuousTab(view: MoleculeRenderView) {

		// restore the wildtype if needed
		posInfo.posConfSpace.wildTypeFragment?.let { setConf(view, it, it.confs.values.first()) }

		// cleanup
		motionInfos.clear()
		stackedMol?.pop()
		stackedMol = null
		renderEffects?.close()
		renderEffects = null
	}

	private sealed class MotionInfo {

		abstract val label: String
		abstract val viewer: MotionViewer?
		abstract fun gui(imgui: Commands)


		class DihedralInfo(
			val view: MoleculeRenderView,
			val motion: ConfLib.ContinuousMotion.DihedralAngle,
			val posInfo: PosInfo,
			val assignmentInfo: Assignments.AssignmentInfo,
			val confConfSpace: ConfSpace.ConfConfSpace,
			val defaultRadiusDegrees: Double
		) : MotionInfo() {

			// make the dihedral angle viewer if this motion is in the conf space
			override var viewer =
				confConfSpace
					.motions
					.filterIsInstance<DihedralAngle.ConfDescription>()
					.find { it.motion === motion }
					?.let { DihedralAngleViewer.make(it, assignmentInfo, view) }

			override val label = run {
				val match = posInfo.pos.findAnchorMatch(confConfSpace.frag)
					?: return@run "Dihedral Angle"
				val atomNames = listOf(motion.a, motion.b, motion.c, motion.d)
					.map { match.resolveName(it) }
				"Dihedral Angle: ${atomNames.joinToString(", ")}"
			}

			override fun gui(imgui: Commands): Unit = imgui.run {

				// show a checkbox to enable/disable the dihedral angle
				checkbox("Included", viewer != null)?.let { isChecked ->
					if (isChecked) {

						// add the motion
						val desc = DihedralAngle.ConfDescription(
							posInfo.pos,
							motion,
							confConfSpace.conf,
							defaultRadiusDegrees
						)
						confConfSpace.motions.add(desc)

						// create the dihedral angle viewer
						viewer = DihedralAngleViewer.make(desc, assignmentInfo, view)

					} else {

						// remove the motion
						confConfSpace.motions.removeIf {
							(it as? DihedralAngle.ConfDescription)?.motion === motion
						}

						// cleanup/reset the viewer
						viewer?.reset()
						viewer = null
					}
				}

				// render the viewer if needed
				viewer?.gui(imgui)
			}
		}
	}
}
