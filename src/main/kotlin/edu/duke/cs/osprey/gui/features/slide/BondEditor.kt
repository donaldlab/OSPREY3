package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.*
import cuchaz.kludge.vulkan.ColorRGBA
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.*
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.molscope.render.HoverEffects
import edu.duke.cs.osprey.molscope.render.MoleculeRenderEffects
import edu.duke.cs.osprey.molscope.render.RenderEffect
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.forcefield.amber.inferBondsAmber
import edu.duke.cs.osprey.gui.prep.BondGuesser
import edu.duke.cs.osprey.gui.prep.covalentRange
import edu.duke.cs.osprey.gui.prep.toTree
import kotlinx.coroutines.runBlocking
import java.util.*


class BondEditor : SlideFeature {

	override val id = FeatureId("edit.bonds")

	private val winState = WindowState()
	private val clickTracker = ClickTracker()
	private var hoverEffects = null as HoverEffects.Writer?
	private val renderEffects = IdentityHashMap<MoleculeRenderView,MoleculeRenderEffects.Writer>()

	private val pMaxDist = Ref.of(3f)

	private var selectedAtom: SelectedAtom? = null

	private var pShowChainBreaks = Ref.of(false)
	private var selectedChain: SelectedChain? = null

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Bonds")) {
			winState.isOpen = true
		}
	}

	private fun Slide.Locked.molViews() = views.mapNotNull { it as? MoleculeRenderView }

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		val views = slide.molViews()

		winState.render(
			onOpen = {

				// add the hover effect
				hoverEffects = slidewin.hoverEffects.writer().apply {
					effect = hoverEffect
				}

				// init the render effects for each molecule
				for (view in views) {
					renderEffects[view] = view.renderEffects.writer()
				}
			},
			whenOpen = {

				// did we click anything?
				if (clickTracker.clicked(slidewin)) {

					selectedAtom = null
					renderEffects.values.forEach { it.clear() }

					// select the atom from the click, if any
					slidewin.mouseTarget?.let { target ->
						(target.view as? MoleculeRenderView)?.let { view ->
							(target.target as? Atom)?.let { atom ->
								selectAtom(view, atom)
							}
						}
					}
				}

				// draw the main window
				window("Bond Editor##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					text("Tools:")
					indent(10f)

					if (button("Clear all bonds")) {
						clearBonds(views)
					}

					if (button("Add bonds automatically")) {
						slidewin.showExceptions {
							guessBonds(views)
						}
					}
					sameLine()
					infoTip("""
						|This tool infers atom connectivity using a molecular mechanics forcefield.
						|After bonds have been added automatically, feel free to use
						|the fine-grained editing tools to add any missing bonds, or
						|remove any extraneous bonds.
					""".trimMargin())

					if (button("Chain Breaks")) {
						pShowChainBreaks.value = true
					}
					sameLine()
					infoTip("""
						|Sometimes the automatic bond tool will add bonds across chain breaks
						|in a polymer, like a protein. This tool allows you to easily remove
						|those incorrectly-added bonds.
					""".trimMargin())

					unindent(10f)

					val selection = selectedAtom

					// show the selected atom
					text("Selected:")
					child("selected", 300f, 30f, true) {
						if (selection != null) {
							text(selection.atom.name)
						} else {
							text("(Click an atom to show bonding options.)")
						}
					}

					// show the nearby atoms
					text("Nearby Atoms:")
					child("nearby", 300f, 300f, true) {
						if (selection != null) {

							columns(2)
							for ((index, info) in selection.nearbyAtoms.withIndex()) {

								// show a checkbox to toggle the bond on/off
								val isCheckHovered = renderBondCheck(imgui, info, index)
								nextColumn()
								renderBondLabel(imgui, info)
								nextColumn()

								// update atom selections
								renderEffects[selection.view]?.set(info.atom, if (isCheckHovered) {
									// highlight the atom when we mouseover the checkbox
									hoverEffect
								} else {
									// otherwise, color by range
									if (info.dist in info.covalentRange) {
										inRangeEffect
									} else {
										outOfRangeEffect
									}
								})
							}
							columns(1)
						}
					}
				}

				// draw the chain breaks window, if needed
				if (pShowChainBreaks.value) {
					window("Bond Editor - Chain Breaks##${slide.name}", pShowChainBreaks, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

						tabBar("polymers") {

							for ((moli, view) in views.withIndex()) {
								val mol = view.currentMol as? Polymer ?: continue

								tabItem("$mol##polymer$moli") {

									for (chain in mol.chains) {

										// pick a default chain if needed
										if (mol.chains.none { it === selectedChain?.chain }) {
											selectedChain = SelectedChain(view, mol, chain)
										}

										if (radioButton(chain.id, selectedChain?.chain === chain)) {
											selectedChain = SelectedChain(view, mol, chain)
										}
									}

									selectedChain?.drawChainBreaks(imgui)
								}
							}
						}
					}
				}
			},
			onClose = {

				// remove the hover effect
				hoverEffects?.close()
				hoverEffects = null

				// clear any leftover selections when the window closes
				renderEffects.values.forEach { it.close() }
				renderEffects.clear()
			}
		)
	}

	override fun contextMenu(contextMenu: ContextMenu, slide: Slide.Locked, slidewin: SlideCommands, target: ViewIndexed) {

		if (!winState.isOpen) {
			return
		}

		// get the view and atom, if any
		val view = target.view as? MoleculeRenderView ?: return
		val atom = target.target as? Atom ?: return

		contextMenu.add {

			// show a button to select the atom
			if (button("Select Atom")) {
				closeCurrentPopup()
				renderEffects.values.forEach { it.clear() }
				selectAtom(view, atom)
			}

			selectedAtom?.let { selection ->
				if (selection.atom !== atom) {

					// get the atom info (reuse a nearby atom if possible)
					val info = selection.nearbyAtoms
						.find { it.atom === atom }
						?: selection.AtomInfo(atom)

					// show a checkbox to toggle the bond on/off
					renderBondCheck(this, info, 0)
					sameLine()
					renderBondLabel(this, info)
				}
			}
		}
	}

	private fun renderBondCheck(imgui: Commands, info: SelectedAtom.AtomInfo, index: Int): Boolean = imgui.run {

		// add a checkbox to toggle the bond
		if (checkbox("${info.atom.name}##atom$index", info.pBonded)) {
			info.updateBond()
		}
		return isItemHovered()
	}

	private fun renderBondLabel(imgui: Commands, info: SelectedAtom.AtomInfo) = imgui.run {

		// show the distance
		textColored(
			if (info.dist in info.covalentRange) {
				ColorRGBA.Int(0, 255, 0)
			} else {
				ColorRGBA.Int(255, 0, 0)
			},
			"${"%.2f".format(info.dist)} A"
		)
		if (isItemHovered()) {
			info.covalentRange.let { range ->
				setTooltip("%.2f A %s in range [%.2f A, %.2f A]".format(
					info.dist,
					if (info.dist in info.covalentRange) {
						"is"
					} else {
						"is not"
					},
					range.start,
					range.endInclusive
				))
			}
		}
	}

	private fun selectAtom(view: MoleculeRenderView, atom: Atom) {

		// make the selection for that atom
		selectedAtom = SelectedAtom(
			view,
			atom,
			pMaxDist.value.toDouble().square()
		).apply {

			// highlight the selected atom
			renderEffects[view]?.set(atom, selectedEffect)
		}
	}

	private fun clearBonds(views: List<MoleculeRenderView>) {
		for (view in views) {
			view.molStack.originalMol.bonds.clear()
			view.moleculeChanged()
		}
	}

	private fun guessBonds(views: List<MoleculeRenderView>) {
		for (view in views) {
			val mol = view.molStack.originalMol
			val bonds = runBlocking { mol.inferBondsAmber() }
			for ((a1, a2) in bonds) {
				mol.bonds.add(a1, a2)
			}
			view.moleculeChanged()
		}
	}
}


private class SelectedAtom(val view: MoleculeRenderView, val atom: Atom, val maxDistSq: Double) {

	val mol get() = view.molStack.originalMol
	val bondGuesser = BondGuesser()

	inner class AtomInfo(val atom: Atom, val dist: Double) {

		constructor(atom: Atom) : this(
			atom,
			this@SelectedAtom.atom.pos.distance(atom.pos)
		)

		val pBonded = Ref.of(mol.bonds.isBonded(this@SelectedAtom.atom, atom))

		val covalentRange = atom.covalentRange(this@SelectedAtom.atom, bondGuesser)

		fun updateBond() {
			val a1 = this@SelectedAtom.atom
			val a2 = this@AtomInfo.atom
			if (pBonded.value) {
				mol.bonds.add(a1, a2)
			} else {
				mol.bonds.remove(a1, a2)
			}
			view.moleculeChanged()
		}
	}

	val nearbyAtoms = mol.atoms.toTree().nearest(atom).asSequence()
		.takeWhile { (_, distSq) -> distSq <= maxDistSq }
		.filter { (nearbyAtom, _) -> atom !== nearbyAtom }
		.toMutableList().apply {
			// add bonded atoms, regardless of distance, so it's easier to break bad bonds
			for (bondedAtom in mol.bonds.bondedAtoms(atom)) {
				add(bondedAtom to atom.pos.distanceSquared(bondedAtom.pos))
			}
		}
		.map { (nearbyAtom, distSq) -> AtomInfo(nearbyAtom, distSq.sqrt()) }
		.toList()
		.sortedBy { it.dist }
}


private class SelectedChain(val view: MoleculeRenderView, val mol: Polymer, val chain: Polymer.Chain) {

	val bondGuesser = BondGuesser()

	// find the connecting bonds for every adjacent pair of residues
	inner class InterBond(
		val res1: Polymer.Residue,
		val res2: Polymer.Residue,
		val bonds: MutableList<Bond>
	) {
		fun Polymer.Residue.label() = "$type $id"
		val label = "${res1.label()} - ${res2.label()}"

		fun removeBonds() {

			// update the molecule
			for (bond in bonds) {
				mol.bonds.remove(bond.atom1, bond.atom2)
			}
			view.moleculeChanged()

			// update this
			bonds.clear()
		}
	}

	inner class Bond(
		val atom1: Atom,
		val atom2: Atom
	) {
		val dist = atom1.pos.distance(atom2.pos)
		val covalentRange = atom1.covalentRange(atom2)
		val good = dist in covalentRange
	}

	private val interBonds: List<InterBond> = ArrayList<InterBond>().apply {
		for (i in 1 until chain.residues.size) {
			val res1 = chain.residues[i - 1]
			val res2 = chain.residues[i]

			val bonds = ArrayList<Bond>()
			for (atom1 in res1.atoms) {
				for (atom2 in mol.bonds.bondedAtoms(atom1)) {
					if (res2.atoms.any { it === atom2 }) {
						bonds.add(Bond(atom1, atom2))
					}
				}
			}

			add(InterBond(res1, res2, bonds))
		}
	}

	fun drawChainBreaks(imgui: Commands) = imgui.run {

		child("chain", 300f, 600f, true) {

			for ((i, interBond) in interBonds.withIndex()) {

				text(interBond.label)
				sameLine()
				if (button("X##$i")) {
					interBond.removeBonds()
				}
				sameLine()
				infoTip("""
					|Remove the bonds between these residues
				""".trimMargin())

				indent(10f)

				// show individual bonds
				if (interBond.bonds.isNotEmpty()) {
					for (bond in interBond.bonds) {

						text("${bond.atom1.name} - ${bond.atom2.name}")
						sameLine()
						textColored(
							if (bond.good) {
								ColorRGBA.Int(0, 255, 0)
							} else {
								ColorRGBA.Int(255, 0, 0)
							},
							"${"%.2f".format(bond.dist)} A"
						)
						if (isItemHovered()) {
							setTooltip("%.2f A %s in range [%.2f A, %.2f A]".format(
								bond.dist,
								if (bond.good) {
									"is"
								} else {
									"is not"
								},
								bond.covalentRange.start,
								bond.covalentRange.endInclusive
							))
						}
					}
				} else {
					text("(no bonds found)")
				}

				unindent(10f)
			}
		}
	}
}


private val hoverEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset, RenderEffect.Flags.Outset),
	200u, 200u, 200u
)
private val selectedEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset, RenderEffect.Flags.Outset),
	255u, 255u, 255u
)
private val inRangeEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset),
	0u, 255u, 0u
)
private val outOfRangeEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset),
	255u, 0u, 0u
)
