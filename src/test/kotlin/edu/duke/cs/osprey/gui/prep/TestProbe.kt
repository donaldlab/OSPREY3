package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.io.OspreyService
import edu.duke.cs.osprey.gui.io.fromPDB
import edu.duke.cs.osprey.gui.io.toPDB
import edu.duke.cs.osprey.gui.io.withService
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.service.services.ClashesRequest
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.shouldBe


class TestProbe : FunSpec({

	context("1cc8") {

		test("protein") {
			withService {

				val mol = Molecule.fromPDB(OspreyGui.getResourceAsString("1cc8.pdb"))
				val results = OspreyService.clashes(ClashesRequest(mol.toPDB()))

				results.groups.getValue("small overlap").run {
					vectors.getValue("orange").size shouldBe 115
					vectors.getValue("red").size shouldBe 32
					vectors.getValue("yellow").size shouldBe 311
					vectors.getValue("yellowtint").size shouldBe 833
				}

				results.groups.getValue("bad overlap").run {
					vectors.getValue("hotpink").size shouldBe 10
				}

				results.groups.getValue("vdw contact").run {
					dots.getValue("sky").size shouldBe 3506
					dots.getValue("green").size shouldBe 2862
					dots.getValue("blue").size shouldBe 4028
					dots.getValue("sea").size shouldBe 2818
				}

				// no H atoms yet, so no Hbonds or salt bridges
			}
		}
	}
})
