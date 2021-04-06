package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.ByteFlags
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.*
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.molscope.render.HoverEffects
import edu.duke.cs.osprey.molscope.render.MoleculeRenderEffects
import edu.duke.cs.osprey.molscope.render.RenderEffect
import edu.duke.cs.osprey.molscope.view.MoleculeRenderStack
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.features.components.ConfLibPicker
import edu.duke.cs.osprey.gui.features.components.DesignPositionEditor
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.prep.*
import java.math.BigInteger
import java.text.NumberFormat
import kotlin.collections.ArrayList


class MutationEditor(val confSpace: ConfSpace) : SlideFeature {

	override val id = FeatureId("edit.mutations")

	private val winState = WindowState()

	private inner class MutEditor(val posInfo: PosInfo) {

		val pos get() = posInfo.pos
		val moltype get() = posInfo.moltype
		val mol= pos.mol

		val winState = WindowState()
			.apply { pOpen.value = true }

		val posEditor = DesignPositionEditor(confSpace, pos)
		val mutationsTabState = Commands.TabState()
		var resetTabSelection = true
		var hoverEffects = null as HoverEffects.Writer?

		var stackedMol: MoleculeRenderStack.StackedMol? = null
		var selectionEffects = null as MoleculeRenderEffects.Writer?

		fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

			val view = slide.views
				.filterIsInstance<MoleculeRenderView>()
				.find { it.molStack.originalMol == mol }
				?: throw Error("can't init design position, molecule has no render view")

			winState.render(
				onOpen = {

					// add the hover effect
					hoverEffects = slidewin.hoverEffects.writer().apply {
						effect = hoverEffect
					}

					// init the pos editor
					posEditor.init(view)

					updateSequenceCounts()
				},
				whenOpen = {

					// draw the window
					window("Edit ${posInfo.pos.name} mutations###${slide.name}_mutEditor_${posInfo.label}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

						// if we need to reset the tab selection, make the flags for the first tab
						fun makeFlags() =
							if (resetTabSelection) {
								resetTabSelection = false
								IntFlags.of(Commands.TabItemFlags.SetSelected)
							} else {
								IntFlags.of(Commands.TabItemFlags.None)
							}

						tabBar("tabs") {
							when (moltype) {
								MoleculeType.Protein -> tabItem("Protein", flags = makeFlags()) {
									posEditor.guiProtein(imgui, slidewin, view)
								}
								// TODO: others?
								else -> Unit
							}
							tabItem("Atoms", flags = makeFlags()) {
								posEditor.guiAtoms(imgui, slidewin, view)
							}
							tabItem(mutationsTabState, "Mutations",
								onActivated = {
									activateMutationsTab()
								},
								whenActive = {
									renderMutationsTab(imgui, view)
								},
								onDeactivated = {
									deactivateMutationsTab(view)
								}
							)
						}
					}
				},
				onClose = {

					// cleanup effects
					hoverEffects?.close()
					hoverEffects = null

					// remove any molecule and effects render overrides
					stackedMol?.pop()
					stackedMol = null
					selectionEffects?.close()
					selectionEffects = null

					// deactivate the mutations tab if it's open
					if (mutationsTabState.wasActive) {
						mutationsTabState.wasActive = false
						deactivateMutationsTab(view)
					}

					// cleanup the pos editor
					posEditor.closed()
					mutEditor = null

					// update the parent window
					this@MutationEditor.run {
						resetInfos()
						updateSequenceCounts()
					}
				}
			)
		}

		private val conflibPicker = ConfLibPicker(confSpace).apply {
			onAdd = { activateMutationsTab() }
		}

		inner class SeqInfo(
			val type: String,
			val label: String,
			val frag: ConfLib.Fragment
		) {
			val pSelected = Ref.of(false)

			// pick an arbitrary conformation from the fragment to show in the GUI
			val conf = frag.confs.values.first()
		}

		private val seqInfos = ArrayList<SeqInfo>()
		private var selectedSeqInfo: SeqInfo? = null

		private fun activateMutationsTab() {

			fun add(type: String, isWildtype: Boolean, frag: ConfLib.Fragment): SeqInfo {

				val label = if (isWildtype) {
					"WildType: $type"
				} else {
					type
				}

				val seqInfo = SeqInfo(type, label, frag).apply {

					// is this fragment selected?
					pSelected.value = pos.confSpace.mutations.contains(type)
				}

				seqInfos.add(seqInfo)
				return seqInfo
			}

			// rebuild the sequence infos
			seqInfos.clear()

			// add the wild type first, if we can
			val wildTypeFrag = pos.confSpace.wildTypeFragment
			val wildTypeInfo = if (wildTypeFrag != null) {
				add(wildTypeFrag.type, true, wildTypeFrag)
			} else {
				null
			}

			// select the wildtype info by default
			selectedSeqInfo = wildTypeInfo

			// add fragments from the libraries
			for (conflib in confSpace.conflibs) {
				conflib.fragments.values
					.filter { pos.isFragmentCompatible(it) }
					.groupBy { it.type }
					.toMutableMap()
					.apply {
						// remove the wildtype if we already have it
						wildTypeFrag?.type?.let { remove(it) }
					}
					.toList()
					.sortedBy { (type, _) -> type }
					.forEach { (type, frag) ->
						add(type, false, frag.first())
					}
			}

			// TODO: collect tags for the fragments? eg, hydrophobic, aromatic
		}

		private fun renderMutationsTab(imgui: Commands, view: MoleculeRenderView) = imgui.run {

			// show the conflib picker
			conflibPicker.render(imgui)

			// show the available mutations
			text("Mutations:")
			sameLine()
			infoTip("""
				|Select the mutations you want to include in your design by clicking the checkboxes.
				|You can temporarily preview a mutation by selecting the radion button next to a mutation.
				|All temporary mutations will be reverted when you're finished with the mutation editor.
			""".trimMargin())
			child("mutations", 300f, 400f, true) {
				if (seqInfos.isNotEmpty()) {
					for (info in seqInfos) {
						if (radioButton("##radio-${info.type}", selectedSeqInfo == info)) {
							selectedSeqInfo = info
							setConf(view, info.frag, info.conf)
						}
						sameLine()
						if (checkbox("${info.label}##check-${info.type}", info.pSelected)) {

							// mark the mutation as included or not in the design position conf space
							if (info.pSelected.value) {
								pos.confSpace.mutations.add(info.type)
							} else {
								pos.confSpace.mutations.remove(info.type)

								// also, remove the confs for this mutation, if any
								pos.confSpace.confs.removeByFragmentType(info.type)
							}

							// update the sequence count
							updateSequenceCounts()
						}
					}
				} else {
					text("(no compatible mutations)")
				}
			}
		}

		private fun deactivateMutationsTab(view: MoleculeRenderView) {

			// remove any molecule and effects render overrides
			stackedMol?.pop()
			stackedMol = null
			selectionEffects?.close()
			selectionEffects = null

			// reset the pos editor to re-draw the render effects
			posEditor.resetInfos()

			selectedSeqInfo = null
		}

		private fun setConf(view: MoleculeRenderView, frag: ConfLib.Fragment, conf: ConfLib.Conf) {

			// make the assignment
			val posAssignment = PosAssignment(posInfo.pos, frag, conf)
			val assignmentInfo = confSpace.assign(posAssignment).assignmentInfos.getValue(posAssignment)
			val mol = assignmentInfo.molInfo.assignedMol

			// override the molecule render view with the new molecule
			stackedMol?.replace(mol)
				?: run {
					stackedMol = view.molStack.push(mol)
				}

			if (selectionEffects == null) {
				selectionEffects = view.renderEffects.writer()
			}

			// add render effects to match the position editor's style
			assignmentInfo.confSwitcher.anchorMatch?.posAnchors?.forEach { anchor ->
				for (atom in anchor.anchorAtoms) {
					selectionEffects?.set(atom, DesignPositionEditor.anchorEffect)
				}
			}
			for (atom in assignmentInfo.confSwitcher.currentAtoms) {
				selectionEffects?.set(atom, DesignPositionEditor.selectedEffect)
			}
		}
	}
	private var mutEditor: MutEditor? = null

	private inner class PosInfo(val pos: DesignPosition, val moltype: MoleculeType, val index: Int) {

		val label = "${pos.name}###pos$index"
		val numSequences = pos.confSpace.mutations.size
	}

	private val DesignPosition.confSpace get() =
		this@MutationEditor.confSpace.positionConfSpaces.getOrMake(this)

	private inner class MolInfo(val molType: MoleculeType, val mol: Molecule) {

		val label = mol.toString()
		val posInfos = ArrayList<PosInfo>()
		var numSequences = BigInteger.ZERO

		fun updateSequenceCount() {
			numSequences = confSpace.countSequences()
		}

		fun makeNewPosition(): PosInfo {

			val positions = confSpace.designPositionsByMol
				.getOrPut(mol) { ArrayList() }

			// choose a default but unique name
			val prefix = "Pos "
			val maxNum = positions
				.mapNotNull {
					it.name
						.takeIf { it.startsWith(prefix) }
						?.substring(prefix.length)
						?.toIntOrNull()
				}
				.max()
				?: 0
			val num = maxNum + 1

			// create the position and add it
			val pos = DesignPosition("$prefix$num", "none", mol)
			positions.add(pos)

			val posInfo = PosInfo(pos, molType, posInfos.size)
			posInfos.add(posInfo)
			return posInfo
		}
	}
	private val molInfos = ArrayList<MolInfo>()

	private fun resetInfos() {

		molInfos.clear()

		for ((moltype, mol) in confSpace.mols) {
			MolInfo(moltype, mol).apply {
				molInfos.add(this)
				confSpace.designPositionsByMol[mol]?.forEachIndexed { index, pos ->
					posInfos.add(PosInfo(pos, moltype, index))
				}
			}
		}
	}

	private var selectedPosInfo = null as PosInfo?

	private val sequenceFormatter = NumberFormat.getIntegerInstance()
		.apply {
			isGroupingUsed = true
		}

	private var numSequences = BigInteger.ZERO

	private fun updateSequenceCounts() {
		molInfos.forEach { it.updateSequenceCount() }
		numSequences = molInfos
			.takeIf { it.isNotEmpty() }
			?.map { it.numSequences }
			?.reduce { a, b -> a*b }
			?: BigInteger.ZERO
	}

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Mutations")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		winState.render(
			onOpen = {

				resetInfos()
				updateSequenceCounts()
			},
			whenOpen = {

				// draw the window
				window("Mutation Editor##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					tabBar("tabs") {

						for (molInfo in molInfos) {
							tabItem(molInfo.label) {

								text("${molInfo.posInfos.size} positions(s)")
								child("positions", 300f, 200f, true) {
									columns(2)
									for (posInfo in molInfo.posInfos) {
										if (selectable(posInfo.label, selectedPosInfo === posInfo)) {
											selectedPosInfo = posInfo
										}
										nextColumn()
										text("${posInfo.numSequences} seqs")
										nextColumn()
									}
									columns(1)
								}

								if (button("Add")) {

									// start the position editor
									mutEditor = MutEditor(molInfo.makeNewPosition())
								}

								sameLine()

								val selectedPosInfo = selectedPosInfo

								enabledIf(selectedPosInfo != null) {
									if (button("Edit")) {

										// sadly the compiler isn't quite smart enough to figure out this can't be null
										// so put in a runtime check to make flow typing work correctly
										selectedPosInfo!!

										// start the position editor
										mutEditor = MutEditor(selectedPosInfo)
									}
								}

								sameLine()

								enabledIf(selectedPosInfo != null) {
									if (button("Remove")) {

										// sadly the compiler isn't quite smart enough to figure out this can't be null
										// so put in a runtime check to make flow typing work correctly
										selectedPosInfo!!

										molInfo.posInfos.remove(selectedPosInfo)
										confSpace.designPositionsByMol[molInfo.mol]?.remove(selectedPosInfo.pos)
										confSpace.positionConfSpaces.remove(selectedPosInfo.pos)
										updateSequenceCounts()
									}
								}

								spacing()

								// show the number of sequences so far
								text("Sequences: ${sequenceFormatter.format(molInfo.numSequences)}")
							}
						}
					}

					// for multiple molecules, show the combined sequence count
					if (molInfos.size > 1) {

						spacing()

						text("Combined Sequences: ${sequenceFormatter.format(numSequences)}")
					}

					// render the position editor, when active
					mutEditor?.gui(imgui, slide, slidewin)
				}
			},
			onClose = {

				// cleanup
				molInfos.clear()
			}
		)
	}
}

private val hoverEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset, RenderEffect.Flags.Outset),
	200u, 200u, 200u
)
