package edu.duke.cs.osprey.gui.prep

import cuchaz.kludge.vulkan.Extent2D
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.WindowCommands
import edu.duke.cs.osprey.molscope.gui.features.slide.CloseSlide
import edu.duke.cs.osprey.molscope.gui.features.slide.MenuRenderSettings
import edu.duke.cs.osprey.molscope.gui.features.slide.CameraTool
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.view.BallAndStick
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.defaultRenderSettings
import edu.duke.cs.osprey.gui.features.slide.*
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.forcefield.amber.partition
import java.util.*
import kotlin.NoSuchElementException


class MoleculePrep(
	win: WindowCommands,
	mols: List<Molecule>
) {

	// partition the molecule into pieces AMBER can understand
	// except, combine all the solvent molecules into one "molecule"
	val partition = mols.partition(combineSolvent = true)

	// TODO: edit the name in the GUI somehow?
	// TODO: persist the name in OMOL somehow?
	var name = mols.firstOrNull()?.name ?: "Prep"

	// include all molecules by default, except solvent
	private val isIncluded = IdentityHashMap<Molecule,Boolean>()
		.apply {
			partition
				.map { (type, mol) ->
					this[mol] = when(type) {
						MoleculeType.Solvent -> false
						else -> true
					}
				}
		}

	fun isIncluded(mol: Molecule) =
		isIncluded[mol] ?: throw NoSuchElementException("mol was not found in this prep")

	fun setIncluded(mol: Molecule, isIncluded: Boolean, slide: Slide.Locked) {

		this.isIncluded[mol] = isIncluded

		// update the views
		val existingView = slide.views.find { it is MoleculeRenderView && it.molStack.originalMol == mol }
		if (isIncluded) {
			if (existingView == null) {
				slide.views.add(BallAndStick(mol))
			}
		} else {
			if (existingView != null) {
				slide.views.remove(existingView)
			}
		}
	}

	fun getIncludedTypedMols(): List<Pair<MoleculeType,Molecule>> =
		partition
			.filter { (_, mol) -> isIncluded[mol] == true }

	fun getIncludedMols(): List<Molecule> =
		getIncludedTypedMols()
			.map { (_, mol) -> mol }

	// make the slide last, since many slide features need to access the prep
	val slide = Slide(name, initialSize = Extent2D(640, 480)).apply {
		lock { s ->

			// make a render view for each molecule in the partition
			for (mol in getIncludedMols()) {
				s.views.add(BallAndStick(mol))
			}
			s.camera.lookAtEverything()

			s.features.menu("File") {
				add(SaveOMOL { getIncludedMols() })
				addSeparator()
				add(ExportPDB { getIncludedMols() })
				add(ExportMol2 { getIncludedMols() to name })
				addSeparator()
				addSpacing(4)
				addSeparator()
				add(CloseSlide(win))
			}
			s.features.menu("View") {
				add(CameraTool())
				add(MenuRenderSettings(defaultRenderSettings))
				add(MoleculeNavigator())
				add(ClashViewer())
			}
			s.features.menu("Prepare") {
				add(FilterTool(this@MoleculePrep))
				add(DuplicateAtomsEditor())
				add(MissingAtomsEditor())
				add(BondEditor())
				add(ProtonationEditor())
				add(MinimizerTool())
			}
		}
		win.addSlide(this)
	}
}
