package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.SharedSpec
import edu.duke.cs.osprey.gui.io.OspreyService
import edu.duke.cs.osprey.gui.io.withService
import edu.duke.cs.osprey.service.services.MoleculeFFInfoRequest
import io.kotlintest.shouldBe


class TestForcefieldInfo : SharedSpec({

	group("benzamidine") {

		test("heavy") {
			withService {

				val results = OspreyService.moleculeFFInfo(MoleculeFFInfoRequest(
					mol2 = OspreyGui.getResourceAsString("benzamidine.gaff2.mol2"),
					ffname = "gaff2"
				))

				results.ffinfo shouldBe OspreyGui.getResourceAsString("benzamidine.frcmod")
			}
		}

		test("with hydrogens") {
			withService {

				val results = OspreyService.moleculeFFInfo(MoleculeFFInfoRequest(
					mol2 = OspreyGui.getResourceAsString("benzamidine.h.gaff2.mol2"),
					ffname = "gaff2"
				))

				results.ffinfo shouldBe OspreyGui.getResourceAsString("benzamidine.h.frcmod")
			}
		}
	}
})
