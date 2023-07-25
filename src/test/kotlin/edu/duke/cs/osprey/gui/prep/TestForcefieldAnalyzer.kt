package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.forcefield.Forcefield
import edu.duke.cs.osprey.gui.io.fromOMOL
import edu.duke.cs.osprey.gui.io.withService
import edu.duke.cs.osprey.molscope.molecule.Molecule
import io.kotest.core.spec.style.FunSpec


class TestForcefieldAnalyzer : FunSpec({

	test("don't crash") {

		val protein = Molecule.fromOMOL(OspreyGui.getResourceAsString("1cc8.protein.omol"))[0]
		val benzamidine = Molecule.fromOMOL(OspreyGui.getResourceAsString("benzamidine.omol"))[0]
		val mols = listOf(protein, benzamidine)

		val ffa = ForcefieldAnalyzer(mols)
		ffa.forcefields.add(Forcefield.Amber96.configure {
			vdwScale = 1.0
		})
		ffa.forcefields.add(Forcefield.EEF1.configure {
			scale = 1.0
		})
		ffa.netCharges[benzamidine]?.netCharge = -1

		val ffp = withService {
			ffa.parameterize()
		}

		for (ff in ffp.forcefields) {
			val energy = ffp.calcEnergy(ff.forcefield)
			println("${ff.forcefield}: $energy")
		}
	}
})
