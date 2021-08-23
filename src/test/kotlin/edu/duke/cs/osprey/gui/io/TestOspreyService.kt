package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.SharedSpec
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.service.services.*
import io.kotlintest.shouldBe
import io.kotlintest.shouldNotBe
import java.nio.file.Paths
import edu.duke.cs.osprey.service.OspreyService as Server


inline fun <T> withService(block: () -> T): T {
	Server.Instance(Paths.get(""), wait = false, useVersionPrefix=true).use {
		UserSettings.serviceProvider = UserSettings.ServiceProvider("localhost", https=false)
		return block()
	}
}

class TestOspreyService : SharedSpec({

	// NOTE: these tests mostly make sure nothing crashes, rather than test for correct output

	test("about") {
		withService {
			OspreyService.about()
		}
	}

	test("missingAtoms") {
		withService {
			val pdb = OspreyGui.getResourceAsString("1cc8.protein.pdb")
			val type = MoleculeType.Protein
			OspreyService.missingAtoms(MissingAtomsRequest(pdb, type.defaultForcefieldNameOrThrow.name))
		}
	}

	group("bonds") {

		test("protein") {
			withService {
				val pdb = OspreyGui.getResourceAsString("1cc8.protein.pdb")
				val type = MoleculeType.Protein
				OspreyService.bonds(BondsRequest(pdb, type.defaultForcefieldNameOrThrow.name))
			}
		}

		test("small mol") {
			withService {
				val pdb = OspreyGui.getResourceAsString("benzamidine.pdb")
				OspreyService.bonds(BondsRequest(pdb, null))
			}
		}
	}

	group("protonation") {

		test("protein") {
			withService {
				val pdb = OspreyGui.getResourceAsString("1cc8.protein.pdb")
				val ffname = MoleculeType.Protein.defaultForcefieldNameOrThrow
				OspreyService.protonation(ProtonationRequest(pdb, ffname.name, null))
			}
		}

		test("small mol") {
			withService {
				val pdb = OspreyGui.getResourceAsString("benzamidine.pdb")
				val ffname = MoleculeType.SmallMolecule.defaultForcefieldNameOrThrow
				OspreyService.protonation(ProtonationRequest(pdb, ffname.name, ffname.atomTypesOrThrow.id))
			}
		}
	}

	group("protonate") {

		test("methylene") {
			withService {

				val c = Atom(Element.Carbon, "C", 0.0, 0.0, 0.0)
				val mol = Molecule(description().name).apply {
					atoms.add(c)
				}

				OspreyService.protonate(ProtonateRequest(
					mol.toMol2(),
					c.name, "c1",
					emptyList(),
					listOf(
						ProtonateRequest.Hydrogen("h1", "hc"),
						ProtonateRequest.Hydrogen("h2", "hc")
					)
				))
			}
		}

		test("diazynediium") {
			withService {

				val n1 = Atom(Element.Nitrogen, "N1",  0.6334, 0.0, 0.0)
				val n2 = Atom(Element.Nitrogen, "N2", -0.6334, 0.0, 0.0)
				val mol = Molecule(description().name).apply {
					atoms.add(n1)
					atoms.add(n2)
					bonds.add(n1, n2)
				}

				OspreyService.protonate(ProtonateRequest(
					mol.toMol2(),
					n1.name, "n1",
					listOf(
						ProtonateRequest.Bond(n2.name, "n1", "T")
					),
					listOf(
						ProtonateRequest.Hydrogen("h1", "hn")
					)
				))
			}
		}
	}

	group("types") {

		test("protein") {
			withService {
				val pdb = OspreyGui.getResourceAsString("1cc8.protein.pdb")
				val ffname = MoleculeType.Protein.defaultForcefieldNameOrThrow
				OspreyService.types(TypesRequest.MoleculeSettings(pdb, ffname.name).toRequest())
			}
		}

		test("small mol") {
			withService {
				val mol2 = OspreyGui.getResourceAsString("benzamidine.h.gaff2.mol2")
				val ffname = MoleculeType.SmallMolecule.defaultForcefieldNameOrThrow
				OspreyService.types(TypesRequest.SmallMoleculeSettings(mol2, ffname.atomTypesOrThrow.id).toRequest())
			}
		}

		test("small mol w/ charges") {
			withService {
				val mol2 = OspreyGui.getResourceAsString("benzamidine.h.gaff2.mol2")
				val ffname = MoleculeType.SmallMolecule.defaultForcefieldNameOrThrow
				val chargeSettings = TypesRequest.ChargeSettings(
					chargeMethod = "bcc",
					netCharge = -1,
					numMinimizationSteps = 0
				)
				OspreyService.types(TypesRequest.SmallMoleculeSettings(mol2, ffname.atomTypesOrThrow.id, chargeSettings).toRequest())
			}
		}
	}

	group("moleculeFFInfo") {

		test("protein") {
			withService {
				val mol2 = OspreyGui.getResourceAsString("1cc8.protein.h.amber.mol2")
				val ffname = MoleculeType.Protein.defaultForcefieldNameOrThrow
				OspreyService.moleculeFFInfo(MoleculeFFInfoRequest(mol2, ffname.name)).ffinfo shouldBe null
			}
		}

		test("benzamidine") {
			withService {
				val mol2 = OspreyGui.getResourceAsString("benzamidine.h.gaff2.mol2")
				val ffname = MoleculeType.SmallMolecule.defaultForcefieldNameOrThrow
				OspreyService.moleculeFFInfo(MoleculeFFInfoRequest(mol2, ffname.name)).ffinfo shouldNotBe null
			}
		}
	}

	group("forcefieldParams") {

		test("protein") {
			withService {
				OspreyService.forcefieldParams(ForcefieldParamsRequest(
					ForcefieldParamsRequest.MolInfo(
						OspreyGui.getResourceAsString("1cc8.protein.h.amber.mol2"),
						MoleculeType.Protein.defaultForcefieldNameOrThrow.name,
						emptyList()
					)
				))
			}
		}

		test("benzamidine") {
			withService {
				OspreyService.forcefieldParams(ForcefieldParamsRequest(
					ForcefieldParamsRequest.MolInfo(
						OspreyGui.getResourceAsString("benzamidine.h.gaff2.mol2"),
						MoleculeType.SmallMolecule.defaultForcefieldNameOrThrow.name,
						listOf(OspreyGui.getResourceAsString("benzamidine.h.frcmod"))
					)
				))
			}
		}

		test("protein and benzamidine") {
			withService {
				OspreyService.forcefieldParams(ForcefieldParamsRequest(
					ForcefieldParamsRequest.MolInfo(
						OspreyGui.getResourceAsString("1cc8.protein.h.amber.mol2"),
						MoleculeType.Protein.defaultForcefieldNameOrThrow.name,
						emptyList()
					),
					ForcefieldParamsRequest.MolInfo(
						OspreyGui.getResourceAsString("benzamidine.h.gaff2.mol2"),
						MoleculeType.SmallMolecule.defaultForcefieldNameOrThrow.name,
						listOf(OspreyGui.getResourceAsString("benzamidine.h.frcmod"))
					)
				))
			}
		}
	}

	group("minimize") {

		test("protein") {
			withService {
				OspreyService.minimize(MinimizeRequest(
					OspreyGui.getResourceAsString("1cc8.protein.top"),
					OspreyGui.getResourceAsString("1cc8.protein.crd"),
					numCycles = 10,
					restraintMask = "@1,2,3,4,5"
				))
			}
		}

		test("benzamidine") {
			withService {
				OspreyService.minimize(MinimizeRequest(
					OspreyGui.getResourceAsString("benzamidine.top"),
					OspreyGui.getResourceAsString("benzamidine.crd"),
					numCycles = 10,
					restraintMask = "@1,2,3,4,5"
				))
			}
		}
	}
})
