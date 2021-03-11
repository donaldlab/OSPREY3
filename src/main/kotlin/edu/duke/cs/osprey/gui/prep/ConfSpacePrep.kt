package edu.duke.cs.osprey.gui.prep

import cuchaz.kludge.vulkan.Extent2D
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.WindowCommands
import edu.duke.cs.osprey.molscope.gui.features.slide.CloseSlide
import edu.duke.cs.osprey.molscope.gui.features.slide.MenuRenderSettings
import edu.duke.cs.osprey.molscope.gui.features.slide.CameraTool
import edu.duke.cs.osprey.molscope.view.BallAndStick
import edu.duke.cs.osprey.gui.defaultRenderSettings
import edu.duke.cs.osprey.gui.features.slide.*


class ConfSpacePrep(
	win: WindowCommands,
	val confSpace: ConfSpace
) {

	// make the slide last, since many slide features need to access the prep
	val slide = Slide(confSpace.name, initialSize = Extent2D(640, 480)).apply {
		lock { s ->

			// make a render view for each molecule
			for ((_, mol) in confSpace.mols) {
				s.views.add(BallAndStick(mol))
			}
			s.camera.lookAtEverything()

			s.features.menu("File") {
				add(SaveConfSpace(confSpace))
				add(SplitConfSpace(confSpace))
				addSeparator()
				add(CompileConfSpace(confSpace))
				addSeparator()
				add(SaveOMOL { confSpace.mols.map { (_, mol) -> mol } })
				add(ExportPDB { confSpace.mols.map { (_, mol) -> mol } })
				add(ExportMol2 { confSpace.mols.map { (_, mol) -> mol } to name })
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
			s.features.menu("Edit") {
				add(ConfSpaceNameEditor(confSpace))
				add(MutationEditor(confSpace))
				add(ConformationEditor(confSpace))
			}
		}
		win.addSlide(this)
	}
}
