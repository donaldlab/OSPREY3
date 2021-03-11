package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.*
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.*
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.molscope.render.Camera
import edu.duke.cs.osprey.molscope.render.MoleculeRenderEffects
import edu.duke.cs.osprey.molscope.render.RenderEffect
import edu.duke.cs.osprey.molscope.tools.identityHashSet
import edu.duke.cs.osprey.molscope.tools.toIdentitySet
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.forcefield.amber.findTypeOrThrow
import edu.duke.cs.osprey.gui.forcefield.amber.findTypes
import org.joml.Vector3d
import java.util.*


class MoleculeNavigator : SlideFeature {

	override val id = FeatureId("view.molecules")

	private val winState = WindowState()

	private val molTypes = IdentityHashMap<Molecule,MoleculeType>()

	private data class Selection(
		var molecule: Molecule? = null,
		var chain: Polymer.Chain? = null,
		var residue: Polymer.Residue? = null,
		var atom: Atom? = null
	) {
		fun clear() {
			molecule = null
			chain = null
			residue = null
			atom = null
		}

		val isSelected get() =
			molecule != null
			|| chain != null
			|| residue != null
			|| atom != null
	}

	private val renderEffects = IdentityHashMap<MoleculeRenderView,MoleculeRenderEffects.Writer>()
	private val viewSelectors = IdentityHashMap<MoleculeRenderView,MoleculeSelector>()
	private val slideHovers = Selection()
	private val guiHovers = Selection()
	private val highlights = Selection()

	private val focusedResidues = IdentityHashMap<Molecule,MutableSet<Polymer.Residue>>()
	private val pShowDist = Ref.of(3.0)

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Molecules")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		val views = slide.views
			.filterIsInstance<MoleculeRenderView>()

		// render the main window
		winState.render(
			onOpen = {

				// init render effects
				for (view in views) {
					renderEffects[view] = view.renderEffects.writer()
				}
				highlights.clear()

				// backup current selectors
				for (view in views) {
					viewSelectors[view] = view.selector
				}
			},
			whenOpen = {

				// look for hovers from the slide
				slideHovers.clear()
				slidewin.mouseTarget?.let {
					val target = it.target
					val view = it.view as? MoleculeRenderView ?: return@let
					val mol = view.currentMol
					when (target) {
						is Atom -> {
							if (slideHovers.atom !== target) {
								slideHovers.atom = target
								slideHovers.molecule = mol

								// find the residue and chain, if possible
								if (mol is Polymer) {
									mol.findChainAndResidue(target)?.let { (chain, res) ->
										slideHovers.chain = chain
										slideHovers.residue = res
									}
								}
							}
						}
					}
				}

				guiHovers.clear()
				window("Molecules##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					for ((moli, view) in views.withIndex()) {

						val mol = view.currentMol
						val type = molTypes.getOrPut(mol) { mol.findTypeOrThrow() }

						withId(moli) {

							// show a label for the molecule that we can hover over
							selectable("$type: $mol", slideHovers.molecule === mol || guiHovers.molecule === mol)
							if (isItemHovered()) {
								guiHovers.molecule = mol
							}

							// show a context menu to center the camera on the molecule
							popupContextItem("centerCamera") {
								if (button("Center Camera")) {
									slidewin.camera.lookAt(mol.atoms, slide)
									slidewin.camera.changed()
									closeCurrentPopup()
								}
							}

							// render a checkbox to temporarily show/hide the molecule
							checkbox("Show", Ref.of(view::isVisible))
							sameLine()
							infoTip("""
								|Uncheck the box to temporarily hide the molecules/atoms.
								|This does not add or remove any molecules/atoms, so hidden
								|molecules/atoms are still there, just not visible.
							""".trimMargin())

							when (type) {
								MoleculeType.Protein -> guiProtein(imgui, slide, slidewin, views, mol as Polymer)
								MoleculeType.SmallMolecule -> guiSmallMolecule(imgui, slide, slidewin, mol)
								else -> Unit
							}
						}

						// let the entries breathe
						spacing()
						spacing()
						separator()
						spacing()
						spacing()
					}

					// show the slide hovers, if any
					text("Hovered:")
					child("hovered", 300f, 80f, border=true) {
						columns(2) {
							column {
								text("Molecule:")
								text("Chain:")
								text("Residue:")
								text("Atom:")
							}
							column {
								text(slideHovers.molecule?.toString() ?: "")
								text(slideHovers.chain?.id ?: "")
								text(slideHovers.residue?.let { "${it.id} ${it.type}" } ?: "")
								text(slideHovers.atom?.let { "${it.name} (${it.element})" } ?: "")
							}
						}
					}
				}

				// update highlights from hovers
				slideHovers.atom?.let { atom ->
					if (highlights.atom !== atom) {

						// highlight just the single atom
						highlights.clear()
						highlights.atom = atom

						renderEffects.values.forEach { it.clear() }
						renderEffects.values
							.find { it.mol === slideHovers.molecule }
							?.let { it[atom] = hoverEffect }
					}
				}
				?: guiHovers.atom?.let { atom ->
					if (highlights.atom !== atom) {

						// highlight just the single atom
						highlights.clear()
						highlights.atom = atom

						renderEffects.values.forEach { it.clear() }
						renderEffects.values
							.find { it.mol === guiHovers.molecule }
							?.let { it[atom] = hoverEffect }
					}
				}
				?: guiHovers.residue?.let { res ->
					if (highlights.residue !== res) {

						// highlight all the atoms in the residue
						highlights.clear()
						highlights.residue = res

						renderEffects.values.forEach { it.clear() }
						renderEffects.values
							.find { it.mol === guiHovers.molecule }
							?.let { it[res.atoms] = hoverEffect }
					}
				}
				?: guiHovers.chain?.let { chain ->
					if (highlights.chain !== chain) {

						// highlight all the atoms in the chain
						highlights.clear()
						highlights.chain = chain

						renderEffects.values.forEach { it.clear() }
						renderEffects.values
							.find { it.mol === guiHovers.molecule }
							?.let { it[chain.residues.flatMap { it.atoms }] = hoverEffect }
					}
				}
				?: guiHovers.molecule?.let { mol ->
					if (highlights.molecule !== mol) {

						// highlight all the atoms in the chain
						highlights.clear()
						highlights.molecule = mol

						renderEffects.values.forEach { it.clear() }
						renderEffects.values
							.find { it.mol === mol }
							?.let { it[MoleculeSelectors.all] = hoverEffect }
					}
				}
				?: run {
					if (highlights.isSelected) {

						// clear any highlights
						highlights.clear()
						renderEffects.values.forEach { it.clear() }
					}
				}
			},
			onClose = {

				// clear highlights
				renderEffects.values.forEach { it.close() }
				renderEffects.clear()

				// restore the molecule selectors
				for (view in views) {
					view.selector = viewSelectors[view] ?: MoleculeSelectors.all
				}
				viewSelectors.clear()

				// other cleanup
				focusedResidues.clear()
			}
		)
	}

	private fun Camera.lookAt(atoms: Collection<Atom>, slide: Slide.Locked) =
		lookAt(atoms.centroid().toFloat(), slide.views)

	private fun Camera.lookAt(atom: Atom, slide: Slide.Locked) =
		lookAt(atom.pos.toFloat(), slide.views)

	private fun Collection<Atom>.centroid() =
		Vector3d().apply {
			forEach { add(it.pos) }
			div(size.toDouble())
		}

	private fun guiProtein(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands, views: List<MoleculeRenderView>, mol: Polymer) = imgui.run {

		fun Polymer.Residue.addToFocus() {
			focusedResidues[mol]?.add(this)
		}

		fun clearFocus() {
			focusedResidues[mol] = identityHashSet()
		}

		fun Polymer.Residue.removeFromFocus() {
			focusedResidues
				.getOrPut(mol) {
					mol.chains
						.flatMap { it.residues }
						.toIdentitySet()
				}
				.remove(this)
		}

		// show all the chains as buttons
		for (chain in mol.chains) {
			treeNode("Residues for chain ${chain.id}") {

				// show all the residues in columns
				val residuesPerCol = 10
				val numCols = chain.residues.size.divideUp(residuesPerCol)
				setNextWindowContentSize(numCols*110f, 0f)
				child("residues", 500f, 260f, true, flags = IntFlags.of(Commands.BeginFlags.HorizontalScrollBar)) {

					columns(numCols, border = true) {
						for (c in 0 until numCols) {
							column {
								for (i in 0 until residuesPerCol) {

									val res = chain.residues.getOrNull(c*residuesPerCol + i) ?: break

									withId(res.id) {

										val label = "${res.id} ${res.type}"

										group {

											// add a checkbox to show/hide the residue
											val isFocused = focusedResidues[mol]?.let { res in it } ?: true
											checkbox("###$label", isFocused)?.let { isChecked ->
												if (isChecked) {
													res.addToFocus()
												} else {
													res.removeFromFocus()
												}
												updateFocusAtoms(views)
											}

											sameLine()

											// show a label for the residue that we can hover over
											selectable(label, slideHovers.residue === res || guiHovers.residue === res)
										}
										if (isItemHovered()) {
											guiHovers.residue = res
											guiHovers.molecule = mol
										}

										// show a context menu for residue-specific ations
										popupContextItem("resPopup") {

											if (button("Center Camera")) {
												slidewin.camera.lookAt(res.atoms, slide)
												slidewin.camera.changed()
												closeCurrentPopup()
											}

											// breathe a little ...
											spacing()
											separator()
											spacing()

											if (button("Show only $label")) {
												clearFocus()
												res.addToFocus()
												updateFocusAtoms(views)
												closeCurrentPopup()
											}

											if (button("Show residues within ${pShowDist.value} Å###showDist")) {
												res
													.findWithin(mol, pShowDist.value)
													.forEach { it.addToFocus() }
												updateFocusAtoms(views)
												closeCurrentPopup()
											}
											sameLine()
											itemWidth(40f) {
												inputDouble("###dist", pShowDist, format = "%.1f")
											}
										}
									}
								}
							}
						}
					}
				}

				if (button("Show all residues")) {
					focusedResidues.remove(mol)
					updateFocusAtoms(views)
				}
				sameLine()
				if (button("Hide all residues")) {
					clearFocus()
					updateFocusAtoms(views)
				}
			}
		}
	}

	private fun guiSmallMolecule(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands, mol: Molecule) = imgui.run {

		treeNode("Atoms") {

			// show all the atoms in columns
			val atomsPerCol = 10
			val numCols = mol.atoms.size.divideUp(atomsPerCol)
			setNextWindowContentSize(numCols*50f, 0f)
			child("atoms", 500f, 200f, true, flags = IntFlags.of(Commands.BeginFlags.HorizontalScrollBar)) {

				columns(numCols, border = true) {
					for (c in 0 until numCols) {
						column {
							for (i in 0 until atomsPerCol) {

								val atom = mol.atoms.getOrNull(c*atomsPerCol + i) ?: break

								withId("$c-$i") {

									// show a label for the atom that we can hover over
									selectable(atom.name, slideHovers.atom === atom || guiHovers.atom === atom)
									if (isItemHovered()) {
										guiHovers.atom = atom
										guiHovers.molecule = mol
									}

									// show a context menu to center the camera on the molecule
									popupContextItem("centerCamera") {
										if (button("Center Camera")) {
											slidewin.camera.lookAt(atom, slide)
											slidewin.camera.changed()
											closeCurrentPopup()
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	override fun contextMenu(contextMenu: ContextMenu, slide: Slide.Locked, slidewin: SlideCommands, target: ViewIndexed) {

		if (!winState.isOpen) {
			return
		}

		// get the atom, if any
		val view = target.view as? MoleculeRenderView ?: return
		val mol = view.currentMol
		val atom = target.target as? Atom ?: return

		contextMenu.add {

			// show details about the atom
			text("Atom: ${atom.name} (${atom.element})")
			indent(10f)

			// show the molecule
			val types = mol.findTypes()
			text("${types.joinToString(",")}: $mol")

			// get the chain and residue, if any
			if (mol is Polymer) {
				mol.findChainAndResidue(atom)?.let { (chain, res) ->
					text("Residue ${chain.id} ${res.id} ${res.type}")
				}
			}

			// show a button to center the camera on the atom
			if (button("Center Camera")) {
				closeCurrentPopup()
				slidewin.camera.lookAt(atom, slide)
				slidewin.camera.changed()
			}

			unindent(10f)
		}
	}

	private fun updateFocusAtoms(views: List<MoleculeRenderView>) {

		for (view in views) {

			val residues = focusedResidues[view.currentMol]
			if (residues == null) {

				// no residue definitions? revert back to default
				view.selector = viewSelectors[view] ?: MoleculeSelectors.all

			} else {

				// otherwise, translate the focus into a selector
				val atoms = residues.flatMap { it.atoms }
				view.selector = { atoms }
			}
		}
	}

	private fun Polymer.Residue.findWithin(mol: Polymer, dist: Double): List<Polymer.Residue> {

		fun dist(a: Polymer.Residue, b: Polymer.Residue) =
			a.atoms.allPairs(b.atoms)
				.map { (a, b) -> a.pos.distance(b.pos) }
				.min()
				?: Double.POSITIVE_INFINITY

		return mol.chains
			.flatMap { it.residues }
			.filter { res -> dist(this, res) <= dist }
	}
}

private val hoverEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Outset),
	200u, 200u, 200u
)
