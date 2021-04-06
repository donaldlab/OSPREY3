package edu.duke.cs.osprey.gui.features.components

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.ByteFlags
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import cuchaz.kludge.tools.divideUp
import edu.duke.cs.osprey.molscope.gui.*
import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.molscope.render.MoleculeRenderEffects
import edu.duke.cs.osprey.molscope.render.RenderEffect
import edu.duke.cs.osprey.molscope.tools.toIdentitySet
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.io.toTomlKey
import edu.duke.cs.osprey.gui.motions.DihedralAngle
import edu.duke.cs.osprey.gui.prep.Anchor
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.gui.prep.DesignPosition
import edu.duke.cs.osprey.gui.prep.Proteins


class DesignPositionEditor(
	val confSpace: ConfSpace,
	val pos: DesignPosition
) {

	val mol = pos.mol

	val dihedralSettings = DihedralAngle.LibrarySettings(
		radiusDegrees = 9.0,
		includeHydroxyls = true,
		includeNonHydroxylHGroups = false
	)

	private val clickTracker = ClickTracker()
	private var autoSelectHydrogens = Ref.of(true)
	private val nameBuf = Commands.TextBuffer(1024)

	private inner class CurrentAtom(val atom: Atom, val label: String) {

		val pSelected = Ref.of(false)
	}

	private fun toggleCurrentAtom(atom: Atom, view: MoleculeRenderView) {

		// is the atom already selected?
		val isSelected = pos.sourceAtoms.any { it === atom }
		if (isSelected) {

			// atom already selected, so deselect it
			pos.sourceAtoms.removeIf { it === atom }

			// toggle the bonded hydrogens if needed
			if (autoSelectHydrogens.value) {
				val hAtoms = mol.bonds.bondedAtoms(atom)
					.filter { it.element == Element.Hydrogen }
					.toIdentitySet()
				pos.sourceAtoms.removeIf { it in hAtoms }
			}

		} else {

			// atom not selected yet, so select it
			pos.sourceAtoms.add(atom)

			// toggle the bonded hydrogens if needed
			if (autoSelectHydrogens.value) {
				mol.bonds.bondedAtoms(atom)
					.filter { it.element == Element.Hydrogen }
					.forEach { h -> pos.sourceAtoms.add(h) }
			}
		}

		resetPosConfSpace()
	}

	private fun removeCurrentAtoms(atoms: Set<Atom>, view: MoleculeRenderView) {

		pos.sourceAtoms.removeIf { it in atoms }

		resetPosConfSpace()
	}

	private inner class AnchorAtom(val atom: Atom, val label: String) {

		var atomClickHandler: ((MoleculeRenderView, Atom) -> Unit)? = null
	}

	private inner class AnchorInfo(
		val anchor: Anchor
	) {

		val anchorAtoms = ArrayList<AnchorAtom>()

		fun findAtomInfo(atom: Atom) =
			anchorAtoms.find { it.atom === atom }
	}

	private inner class AnchorGroupInfo(
		val anchors: MutableList<Anchor>
	) {

		val anchorInfos = ArrayList<AnchorInfo>()

		fun replaceAnchor(old: Anchor, new: Anchor) {
			anchors[anchors.indexOf(old)] = new
		}
	}

	private val currentAtoms = ArrayList<CurrentAtom>()
	private val anchorGroupInfos = ArrayList<AnchorGroupInfo>()
	private var renderEffects = null as MoleculeRenderEffects.Writer?

	fun init(view: MoleculeRenderView) {

		renderEffects = view.renderEffects.writer()

		resetInfos()
	}

	fun closed() {

		renderEffects?.close()
		renderEffects = null

		anchorGroupInfos.clear()
		currentAtoms.clear()
	}

	fun Atom.label(): String {

		// use the atom name as the label
		var label = name

		// if there's a residue, append that too
		if (mol is Polymer) {
			val res = mol.findResidue(this)
			if (res != null) {
				label += " @ $res"
			}
		}

		return label
	}

	private fun resetPosConfSpace() {

		// delete the old pos conf space and make a new one
		confSpace.positionConfSpaces.remove(pos)
		confSpace.positionConfSpaces.getOrMake(pos)

		// start with the wild-type "mutation" by default
		confSpace.addMutations(pos, pos.type)

		// add the libary confs
		confSpace.addConformationsFromLibraries(pos, pos.type)

		// add the wild-type frag/conf
		try {
			confSpace.addWildTypeConformation(pos)
		} catch (ex: DesignPosition.IllegalAnchorsException) {
			// ignore this error here
			// it means there's something wrong with the anchor configuration in the conflib,
			// but we don't want to take down the whole GUI with a crash
		}

		// add default dihedral angles from the conformation library
		for (type in confSpace.getMutations(pos)) {
			for (confConfSpace in confSpace.getConformations(pos, type)) {
				for (motion in DihedralAngle.ConfDescription.makeFromLibrary(pos, confConfSpace.frag, confConfSpace.conf, dihedralSettings)) {
					confConfSpace.motions.add(motion)
				}
			}
		}

		resetInfos()
	}

	fun resetInfos() {

		renderEffects?.clear()
		anchorGroupInfos.clear()
		currentAtoms.clear()

		// make infos for the anchor atoms
		for (anchorGroup in pos.anchorGroups) {
			anchorGroupInfos.add(AnchorGroupInfo(anchorGroup).apply {
				for (anchor in anchorGroup) {
					anchorInfos.add(AnchorInfo(anchor).apply {
						for (atom in anchor.anchorAtoms) {
							val info = AnchorAtom(atom, atom.label())
							anchorAtoms.add(info)
							renderEffects?.set(atom, anchorEffect)
						}
					})
				}
			})
		}

		// make infos for the current atoms
		for (atom in pos.sourceAtoms) {
			val info = CurrentAtom(atom, atom.label())
			currentAtoms.add(info)
			renderEffects?.set(atom, selectedEffect)
		}
	}

	private fun <R> hoveredAtom(slidewin: SlideCommands, block: (MoleculeRenderView, Atom) -> R): R? {

		// select the atom from the click, if any
		slidewin.mouseTarget?.let { target ->
			(target.view as? MoleculeRenderView)?.let { view ->

				// make sure we're picking from the same molecule
				if (view.molStack.originalMol == mol) {

					(target.target as? Atom)?.let { atom ->

						// found an elligible atom!
						return block(view, atom)
					}
				}
			}
		}

		// no elligible atom hovered
		return null
	}

	private var selectedRes: Polymer.Residue? =
		// try to guess the selected residue from the design position
		(pos.mol as? Polymer)?.let { mol ->
			pos.sourceAtoms
				.mapNotNull { mol.findResidue(it) }
				.toIdentitySet()
				.takeIf { it.size == 1 }
				?.firstOrNull()
		}

	private var selectedChain: Polymer.Chain? =
		// try to guess the selected chain from the selected residue
		selectedRes?.let { res ->
			val mol = pos.mol as Polymer
			mol.chains.find { chain -> chain.residues.any { it === res } }
		}

	fun guiProtein(imgui: Commands, slidewin: SlideCommands, view: MoleculeRenderView) = imgui.run {

		if (mol !is Polymer) return

		// edit the name
		nameBuf.text = pos.name
		if (inputText("Name", nameBuf)) {
			pos.name = nameBuf.text
		}

		// select an initial chain if needed
		if (selectedChain == null) {
			selectedChain = mol.chains.first()
		}

		// did we pointing at an atom?
		val hoveredRes = hoveredAtom(slidewin) { view, atom ->

			// was it a residue?
			val res = mol.findResidue(atom)
			if (res != null) {

				// did we click on it?
				if (clickTracker.clicked(slidewin)) {
					selectResidue(view, res)
				}
			}

			res
		}

		// show all the chains as radios
		text("Chains:")
		columns(5)
		for (chain in mol.chains) {
			if (radioButton(chain.id, selectedChain == chain)) {
				selectedChain = chain
			}
			nextColumn()
		}
		columns(1)

		// we should definitely have a chain by now
		val selectedChain = selectedChain ?: return

		// show all the residues as radios
		val residuesPerCol = 20
		val numCols = selectedChain.residues.size.divideUp(residuesPerCol)
		text("Residues:")
		setNextWindowContentSize(numCols*100f, 0f)
		child("residues", 500f, 490f, true, flags = IntFlags.of(Commands.BeginFlags.HorizontalScrollBar)) {

			columns(numCols, border = true) {
				for (c in 0 until numCols) {
					column {
						for (i in 0 until residuesPerCol) {

							val res = selectedChain.residues.getOrNull(c*residuesPerCol + i) ?: break

							if (radioButton("${res.id} ${res.type}", res == selectedRes || res == hoveredRes)) {
								slidewin.showExceptions {
									selectResidue(view, res)
								}
							}
						}
					}
				}
			}
		}
	}

	private fun selectResidue(view: MoleculeRenderView, res: Polymer.Residue) {

		// save the selection
		selectedRes = res

		// update the design position
		Proteins.setDesignPosition(pos, res)
		resetPosConfSpace()
	}

	private var atomClickHandler: ((MoleculeRenderView, Atom) -> Unit)? = null

	private fun atomClickSidechain(view: MoleculeRenderView, atom: Atom) {
		toggleCurrentAtom(atom, view)
	}

	fun guiAtoms(imgui: Commands, slidewin: SlideCommands, view: MoleculeRenderView) = imgui.run {

		var posChanged = false

		// default to the sidechain atom click handler
		if (atomClickHandler == null) {
			atomClickHandler = ::atomClickSidechain
		}

		// did we click on an atom?
		hoveredAtom(slidewin) { view, atom ->
			if (clickTracker.clicked(slidewin)) {

				// invoke the atom click handler
				atomClickHandler?.invoke(view, atom)
			}
		}

		// edit the name
		nameBuf.text = pos.name
		if (inputText("Name", nameBuf)) {
			pos.name = nameBuf.text
		}

		spacing()

		// show current atoms (but call them "sidechain" atoms I guess?)
		if (radioButton("Sidechain Atoms: ${currentAtoms.size}", atomClickHandler == ::atomClickSidechain)) {
			atomClickHandler = ::atomClickSidechain
		}
		indent(20f)
		checkbox("Automatically select hydrogens", autoSelectHydrogens)
		child("currentAtoms", 300f, 200f, true) {
			for (info in currentAtoms) {
				selectable(info.label, info.pSelected)
			}
		}

		enabledIf(currentAtoms.any { it.pSelected.value }) {
			if (button("Remove")) {
				val selectedAtoms = currentAtoms
					.filter { it.pSelected.value }
					.map { it.atom }
					.toIdentitySet()
				removeCurrentAtoms(selectedAtoms, view)
			}
		}

		unindent(20f)
		spacing()

		// TODO: cleanup anchor rendering code? the indentation is getting a little ridiculous...

		// show anchors
		text("Anchors:")
		sameLine()
		infoTip("""
				|Anchor atoms are mainchain atoms that help align the sidechain to the mainchain.
			""".trimMargin())
		indent(20f)

		for ((anchorsi, anchorGroupInfo) in anchorGroupInfos.withIndex()) {

			// breathe a little
			if (anchorsi > 0) {
				spacing()
				separator()
				spacing()
			}

			text("Anchor Group:")
			sameLine()
			withId(anchorsi) {

				// button to delete anchor group
				if (button("x")) {
					pos.anchorGroups.remove(anchorGroupInfo.anchors)
					posChanged = true
				}

				for ((anchori, anchorInfo) in anchorGroupInfo.anchorInfos.withIndex()) {
					indent(20f)
					withId(anchori) {

						fun anchorDeleteButton() {
							if (button("x")) {
								anchorGroupInfo.anchors.remove(anchorInfo.anchor)
								posChanged = true
							}
						}

						fun anchorRadioButton(name: String, atom: Atom, anchorUpdater: (Atom) -> Anchor) {
							val atomInfo = anchorInfo.findAtomInfo(atom)
							if (atomInfo == null) {
								text("$name: (error)")
							} else {
								// TODO: BUGBUG: having the same atom in multiple spots causes all the radios to light up
								if (radioButton("$name: ${atomInfo.label}", atomClickHandler === atomInfo.atomClickHandler)) {

									/* This is a little tricky, so let me explain:

										When the radio is clicked, we make a function that:
											1. creates new anchor with the selected atom by calling anchorUpdater()
											2. replaces the old anchor with the new anchor in the anchor group
											3. rebuilds all the GUI shadow info for the design position by calling resetInfos()
											4. resets the atoms tab atom click handler

										Then we assign that function to the atom info for the anchor atom so we can find it again later.
										And we also assign it to the atoms tab atom click handler, so it's activated.

										Conveniently, resetInfos() will drop this function from the heap when
										it re-creates all the AnchorInfos instances, so we won't have stale handlers hanging around
									*/
									atomInfo.atomClickHandler = { view, clickedAtom ->
										anchorGroupInfo.replaceAnchor(
											anchorInfo.anchor,
											anchorUpdater(clickedAtom)
										)
										/* NOTE:
											Could set posChanged = true here, but it won't work the way you'd think.
											By the time atomClickHandler gets called, we'll be on a different stack frame,
											and a different instance of posChanged.
											So just call resetInfos() directly. Which is ok here, since
											atomClickHandler gets called outside of the info render loops.
										 */
										resetPosConfSpace()
										atomClickHandler = null
									}
									atomClickHandler = atomInfo.atomClickHandler
								}
							}
						}

						when (anchorInfo.anchor) {
							is Anchor.Single -> {
								text("Single Anchor")
								sameLine()
								infoTip("""
										|Single anchors allow the sidechain to bond to the mainchain at one atom.
										|Three atoms are required: a, b, and c.
										|Only the a atom is allowed to bond to the sidechain.
										|First, a positions are exactly superposed.
										|Then, a->b vectors are made parallel.
										|Then, a,b,c planes are made parallel.
									""".trimMargin())
								sameLine()
								anchorDeleteButton()

								anchorRadioButton("a", anchorInfo.anchor.a) { pickedAtom ->
									anchorInfo.anchor.copy(a = pickedAtom)
								}
								anchorRadioButton("b", anchorInfo.anchor.b) { pickedAtom ->
									anchorInfo.anchor.copy(b = pickedAtom)
								}
								anchorRadioButton("c", anchorInfo.anchor.c) { pickedAtom ->
									anchorInfo.anchor.copy(c = pickedAtom)
								}
							}
							is Anchor.Double -> {
								text("Double Anchor")
								sameLine()
								infoTip("""
										|Double anchors allow the sidechain to bond to the mainchain at two atoms.
										|Four atoms are required: a, b, c, and d.
										|Only the a and b atoms are allowed to bond to the sidechain.
										|TODO: describe alignment
									""".trimMargin())
								sameLine()
								anchorDeleteButton()

								anchorRadioButton("a", anchorInfo.anchor.a) { pickedAtom ->
									anchorInfo.anchor.copy(a = pickedAtom)
								}
								anchorRadioButton("b", anchorInfo.anchor.b) { pickedAtom ->
									anchorInfo.anchor.copy(b = pickedAtom)
								}
								anchorRadioButton("c", anchorInfo.anchor.c) { pickedAtom ->
									anchorInfo.anchor.copy(c = pickedAtom)
								}
								anchorRadioButton("d", anchorInfo.anchor.d) { pickedAtom ->
									anchorInfo.anchor.copy(d = pickedAtom)
								}
							}
						}
					}

					unindent(20f)
				}

				if (button("Add Single Anchor")) {
					slidewin.showExceptions {
						anchorGroupInfo.anchors.add(pos.anchorSingle(
							// pick dummy atoms for now
							a = mol.atoms[0],
							b = mol.atoms[1],
							c = mol.atoms[2]
						))
						posChanged = true
					}
				}

				sameLine()

				if (button("Add Double Anchor")) {
					slidewin.showExceptions {
						anchorGroupInfo.anchors.add(pos.anchorDouble(
							// pick dummy atoms for now
							a = mol.atoms[0],
							b = mol.atoms[1],
							c = mol.atoms[2],
							d = mol.atoms[3]
						))
						posChanged = true
					}
				}
			}
		}

		unindent(20f)

		spacing()

		if (button("Add Anchor Group")) {
			pos.anchorGroups.add(ArrayList())
			posChanged = true
		}

		// finally, handle any pending updates after we changed the position
		if (posChanged) {
			resetPosConfSpace()
		}
	}

	companion object {

		val selectedEffect = RenderEffect(
			ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset, RenderEffect.Flags.Outset),
			255u, 255u, 255u
		)
		val anchorEffect = RenderEffect(
			ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset, RenderEffect.Flags.Outset),
			0u, 255u, 0u
		)
	}
}
