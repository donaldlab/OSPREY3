package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.gui.io.toMol2
import edu.duke.cs.osprey.gui.io.withService
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule
import io.kotest.assertions.withClue
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.shouldBe


class TestProtonation : FunSpec({

	fun Molecule.findProtonation(atom: Atom, filter: (Protonation) -> Boolean): Protonation {
		val protonations = protonations(atom)
			.filter(filter)
		when (protonations.size) {
			0 -> throw NoSuchElementException("no protonations matched")
			1 -> return protonations.first()
			else -> throw NoSuchElementException("multiple protonations:" + protonations.joinToString("\n"))
		}
	}

	fun Molecule.shouldHaveUniqueHydrogens() {
		val names = atoms
			.filter { it.element == Element.Hydrogen }
			.map { it.name }
		withClue(names) {
			names.toSet().size shouldBe names.size
		}
	}

	fun Molecule.numH(atom: Atom) = bonds.bondedAtoms(atom).count { it.element == Element.Hydrogen }

	context("carbon") {

		test("methylene") {
			withService {

				val c = Atom(Element.Carbon, "C", 0.0, 0.0, 0.0)
				val mol = Molecule("methylene").apply {
					atoms.add(c)
				}

				val protonation = mol.findProtonation(c) { it.numH == 2 }
				mol.protonate(c, protonation)
				mol.numH(c) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check the H-C-H angle?
			}
		}

		test("methyl cation") {
			withService {

				val c = Atom(Element.Carbon, "C", 0.0, 0.0, 0.0)
				val mol = Molecule("methyl cation").apply {
					atoms.add(c)
				}

				val protonation = mol.findProtonation(c) { it.numH == 3 }
				mol.protonate(c, protonation)
				mol.numH(c) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for trigonal planarity?
			}
		}

		test("methane") {
			withService {

				val c = Atom(Element.Carbon, "C", 0.0, 0.0, 0.0)
				val mol = Molecule("methane").apply {
					atoms.add(c)
				}

				val protonation = mol.findProtonation(c) { it.numH == 4 }
				mol.protonate(c, protonation)
				mol.numH(c) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for tetrahedral geometry?
			}
		}

		test("ethyne") {
			withService {

				val c1 = Atom(Element.Carbon, "C1", 0.0, 0.0, 0.0)
				val c2 = Atom(Element.Carbon, "C2", 0.0, 0.0, 1.1983)
				val mol = Molecule("ethyne").apply {
					atoms.add(c1)
					atoms.add(c2)
					bonds.add(c1, c2)
				}

				val protonation = mol.findProtonation(c1) { it.numH == 1 }
				mol.protonate(c1, protonation)
				mol.numH(c1) shouldBe protonation.numH
				mol.protonate(c2, protonation)
				mol.numH(c2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for linearity?
			}
		}

		test("ethene") {
			withService {

				val c1 = Atom(Element.Carbon, "C1", 0.0, 0.0, 0.0)
				val c2 = Atom(Element.Carbon, "C2", 0.0, 0.0, 1.3343)
				val mol = Molecule("ethene").apply {
					atoms.add(c1)
					atoms.add(c2)
					bonds.add(c1, c2)
				}

				val protonation = mol.findProtonation(c1) { it.numH == 2 }
				mol.protonate(c1, protonation)
				mol.numH(c1) shouldBe protonation.numH
				mol.protonate(c2, protonation)
				mol.numH(c2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for planarity?
			}
		}

		test("ethane") {
			withService {

				val c1 = Atom(Element.Carbon, "C1", 0.0, 0.0, 0.0)
				val c2 = Atom(Element.Carbon, "C2", 0.0, 0.0, 1.5375)
				val mol = Molecule("ethane").apply {
					atoms.add(c1)
					atoms.add(c2)
					bonds.add(c1, c2)
				}

				val protonation = mol.findProtonation(c1) { it.numH == 3 }
				mol.protonate(c1, protonation)
				mol.numH(c1) shouldBe protonation.numH
				mol.protonate(c2, protonation)
				mol.numH(c2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for tetrahedral geometry?
				// TODO: check for staggered conformation?
			}
		}

		test("benzene") {
			withService {

				val c1 = Atom(Element.Carbon, "C1", -2.553, 3.433, 3.188)
				val c2 = Atom(Element.Carbon, "C2", -1.337, 3.205, 2.508)
				val c3 = Atom(Element.Carbon, "C3", -0.791, 4.128, 1.629)
				val c4 = Atom(Element.Carbon, "C4", -1.550, 5.259, 1.337)
				val c5 = Atom(Element.Carbon, "C5", -2.797, 5.495, 1.902)
				val c6 = Atom(Element.Carbon, "C6", -3.296, 4.579, 2.847)
				val mol = Molecule("benzene").apply {
					atoms.add(c1)
					atoms.add(c2)
					bonds.add(c1, c2)
					atoms.add(c3)
					bonds.add(c2, c3)
					atoms.add(c4)
					bonds.add(c3, c4)
					atoms.add(c5)
					bonds.add(c4, c5)
					atoms.add(c6)
					bonds.add(c5, c6)
					bonds.add(c6, c1)
				}

				val protonation = mol.findProtonation(c1) { it.numH == 1 }
				mol.protonate(c1, protonation)
				mol.numH(c1) shouldBe protonation.numH
				mol.protonate(c2, protonation)
				mol.numH(c2) shouldBe protonation.numH
				mol.protonate(c3, protonation)
				mol.numH(c3) shouldBe protonation.numH
				mol.protonate(c4, protonation)
				mol.numH(c4) shouldBe protonation.numH
				mol.protonate(c5, protonation)
				mol.numH(c5) shouldBe protonation.numH
				mol.protonate(c6, protonation)
				mol.numH(c6) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for planar geometry?
			}
		}

		test("cyclohexane") {
			withService {

				val c1 = Atom(Element.Carbon, "C1", 20.5876, -5.5100, -0.5150)
				val c2 = Atom(Element.Carbon, "C2", 21.3226, -4.2453, -0.0471)
				val c3 = Atom(Element.Carbon, "C3", 20.5992, -2.9976, -0.5750)
				val c4 = Atom(Element.Carbon, "C4", 19.1520, -2.9786, -0.0612)
				val c5 = Atom(Element.Carbon, "C5", 18.4174, -4.2435, -0.5289)
				val c6 = Atom(Element.Carbon, "C6", 19.1407, -5.4908,  0.0000)
				val mol = Molecule("cyclohexane").apply {
					atoms.add(c1)
					atoms.add(c2)
					bonds.add(c1, c2)
					atoms.add(c3)
					bonds.add(c2, c3)
					atoms.add(c4)
					bonds.add(c3, c4)
					atoms.add(c5)
					bonds.add(c4, c5)
					atoms.add(c6)
					bonds.add(c5, c6)
					bonds.add(c6, c1)
				}

				val protonation = mol.findProtonation(c1) { it.numH == 1 }
				mol.protonate(c1, protonation)
				mol.numH(c1) shouldBe protonation.numH
				mol.protonate(c2, protonation)
				mol.numH(c2) shouldBe protonation.numH
				mol.protonate(c3, protonation)
				mol.numH(c3) shouldBe protonation.numH
				mol.protonate(c4, protonation)
				mol.numH(c4) shouldBe protonation.numH
				mol.protonate(c5, protonation)
				mol.numH(c5) shouldBe protonation.numH
				mol.protonate(c6, protonation)
				mol.numH(c6) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check tetrahedral geometry?
			}
		}

		test("isobutane") {
			withService {

				val c  = Atom(Element.Carbon, "C",   1.4204, -0.2130, -0.3689)
				val c1 = Atom(Element.Carbon, "C1", -0.6743,  1.2682, -0.3689)
				val c2 = Atom(Element.Carbon, "C2", -0.6743, -0.9536,  0.9139)
				val c3 = Atom(Element.Carbon, "C3", -0.1300, -0.1839, -0.3185)
				val mol = Molecule("isobutane").apply {
					atoms.add(c)
					atoms.add(c1)
					bonds.add(c, c1)
					atoms.add(c2)
					bonds.add(c, c2)
					atoms.add(c3)
					bonds.add(c, c3)
				}

				val protonation = mol.findProtonation(c) { it.numH == 1 }
				mol.protonate(c, protonation)
				mol.numH(c) shouldBe protonation.numH

				val protonationMethyl = mol.findProtonation(c1) { it.numH == 3 }
				mol.protonate(c1, protonationMethyl)
				mol.numH(c1) shouldBe protonationMethyl.numH
				mol.protonate(c2, protonationMethyl)
				mol.numH(c2) shouldBe protonationMethyl.numH
				mol.protonate(c3, protonationMethyl)
				mol.numH(c3) shouldBe protonationMethyl.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check tetrahedral geometry?
			}
		}
	}

	context("nitrogen") {

		test("azanide anion") {
			withService {

				val n = Atom(Element.Nitrogen, "N", 0.0, 0.0, 0.0)
				val mol = Molecule("nitrogen").apply {
					atoms.add(n)
				}

				val protonation = mol.findProtonation(n) { it.numH == 2 }
				mol.protonate(n, protonation)
				mol.numH(n) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check H-N-H angle?
			}
		}

		test("ammonia") {
			withService {

				val n = Atom(Element.Nitrogen, "N", 0.0, 0.0, 0.0)
				val mol = Molecule("ammonia").apply {
					atoms.add(n)
				}

				val protonation = mol.findProtonation(n) { it.numH == 3 }
				mol.protonate(n, protonation)
				mol.numH(n) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check trigonal pyramidal geometry?
			}
		}

		test("ammonium cation") {
			withService {

				val n = Atom(Element.Nitrogen, "N", 0.0, 0.0, 0.0)
				val mol = Molecule("ammonium cation").apply {
					atoms.add(n)
				}

				val protonation = mol.findProtonation(n) { it.numH == 4 }
				mol.protonate(n, protonation)
				mol.numH(n) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check tetrahedral geometry?
			}
		}

		test("diazynediium") {
			withService {

				val n1 = Atom(Element.Nitrogen, "N1",  0.5510, 0.0, 0.0)
				val n2 = Atom(Element.Nitrogen, "N2", -0.5510, 0.0, 0.0)
				val mol = Molecule("diazynediium").apply {
					atoms.add(n1)
					atoms.add(n2)
					bonds.add(n1, n2)
				}

				val protonation = mol.findProtonation(n1) { it.numH == 1 && it.hybridization == Hybridization.Sp }
				mol.protonate(n1, protonation)
				mol.numH(n1) shouldBe protonation.numH
				mol.protonate(n2, protonation)
				mol.numH(n2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check linearity?
			}
		}

		test("diazene") {
			withService {

				val n1 = Atom(Element.Nitrogen, "N1",  0.6100, 0.0, 0.0)
				val n2 = Atom(Element.Nitrogen, "N2", -0.6100, 0.0, 0.0)
				val mol = Molecule("diazene").apply {
					atoms.add(n1)
					atoms.add(n2)
					bonds.add(n1, n2)
				}

				val protonation = mol.findProtonation(n1) { it.numH == 1 && it.hybridization == Hybridization.Sp2 }
				mol.protonate(n1, protonation)
				mol.numH(n1) shouldBe protonation.numH
				mol.protonate(n2, protonation)
				mol.numH(n2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check planar geometry?
				// TODO: what about stereoisomers?
			}
		}

		test("formamide") {
			withService {

				val n = Atom(Element.Nitrogen, "N", -0.3618, 0.5603, 0.0000)
				val c = Atom(Element.Carbon, "C", 0.2763, -0.5443, -0.0000)
				val o = Atom(Element.Oxygen, "O", 1.4963, -0.5453, 0.0000)
				val mol = Molecule("formamide").apply {
					atoms.add(n)
					atoms.add(c)
					bonds.add(n, c)
					atoms.add(o)
					bonds.add(c, o)
				}

				val protonationN = mol.findProtonation(n) { it.numH == 2 && it.hybridization == Hybridization.Sp2 }
				mol.protonate(n, protonationN)
				mol.numH(n) shouldBe protonationN.numH

				val protonationC = mol.findProtonation(c) { it.numH == 1 }
				mol.protonate(c, protonationC)
				mol.numH(c) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check planar geometry
			}
		}

		test("hydrazine") {
			withService {

				val n1 = Atom(Element.Nitrogen, "N1",  0.6334, 0.0, 0.0)
				val n2 = Atom(Element.Nitrogen, "N2", -0.6334, 0.0, 0.0)
				val mol = Molecule("hydrazine").apply {
					atoms.add(n1)
					atoms.add(n2)
					bonds.add(n1, n2)
				}

				val protonation = mol.findProtonation(n1) { it.numH == 2 && it.hybridization == Hybridization.Sp3 }
				mol.protonate(n1, protonation)
				mol.numH(n1) shouldBe protonation.numH
				mol.protonate(n2, protonation)
				mol.numH(n2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check trigonal pyramidal geometry?
			}
		}

		test("diazanediium") {
			withService {

				val n1 = Atom(Element.Nitrogen, "N1",  0.6334, 0.0, 0.0)
				val n2 = Atom(Element.Nitrogen, "N2", -0.6334, 0.0, 0.0)
				val mol = Molecule("diazanediium").apply {
					atoms.add(n1)
					atoms.add(n2)
					bonds.add(n1, n2)
				}

				val protonation = mol.findProtonation(n1) { it.numH == 3 }
				mol.protonate(n1, protonation)
				mol.numH(n1) shouldBe protonation.numH
				mol.protonate(n2, protonation)
				mol.numH(n2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check tetrahedral geometry?
				// TODO: should we expect a staggered conformation?

				println(mol.toMol2())
			}
		}

		test("N-methylformamide") {
			withService {

				val n = Atom(Element.Nitrogen, "N", -0.3662, -0.3144, -0.0000)
				val c = Atom(Element.Carbon, "C", 1.0230, -0.3844, 0.0000)
				val co = Atom(Element.Carbon, "CO", -1.0788, 0.7609, -0.0000)
				val o = Atom(Element.Oxygen, "O", -2.2968, 0.6774, -0.0000)
				val mol = Molecule("N-methylformamide").apply {
					atoms.add(n)
					atoms.add(c)
					bonds.add(n, c)
					atoms.add(co)
					bonds.add(n, co)
					atoms.add(o)
					bonds.add(co, o)
				}

				val protonationN = mol.findProtonation(n) { it.numH == 1 && it.hybridization == Hybridization.Sp2 }
				mol.protonate(n, protonationN)
				mol.numH(n) shouldBe protonationN.numH

				val protonationC = mol.findProtonation(c) { it.numH == 3 }
				mol.protonate(c, protonationC)
				mol.numH(c) shouldBe protonationC.numH

				val protonationCO = mol.findProtonation(co) { it.numH == 1 }
				mol.protonate(co, protonationCO)
				mol.numH(co) shouldBe protonationCO.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check planar geometry
			}
		}

		test("dimethlamine") {
			withService {

				val n = Atom(Element.Nitrogen, "N", -0.3924, -0.5549, -0.0984)
				val c1 = Atom(Element.Carbon, "C1", 1.0787, -0.6370, -0.0828)
				val c2 = Atom(Element.Carbon, "C2", -0.9601, 0.8047, -0.0828)
				val mol = Molecule("dimethlamine").apply {
					atoms.add(n)
					atoms.add(c1)
					bonds.add(n, c1)
					atoms.add(c2)
					bonds.add(n, c2)
				}

				val protonationN = mol.findProtonation(n) { it.numH == 1 && it.hybridization == Hybridization.Sp3 }
				mol.protonate(n, protonationN)
				mol.numH(n) shouldBe protonationN.numH

				val protonationC = mol.findProtonation(c1) { it.numH == 3 }
				mol.protonate(c1, protonationC)
				mol.numH(c1) shouldBe protonationC.numH
				mol.protonate(c2, protonationC)
				mol.numH(c2) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for trigonal pyrimidal geometry?
			}
		}

		test("dimethylammonium cation") {
			withService {

				val n = Atom(Element.Nitrogen, "N", 1.3000, 0.8151, -0.0008)
				val c1 = Atom(Element.Carbon, "C1", 0.0997, -0.0335, 0.0000)
				val c2 = Atom(Element.Carbon, "C2", 2.5006, -0.0331, 0.0002)
				val mol = Molecule("dimethylammonium cation").apply {
					atoms.add(n)
					atoms.add(c1)
					bonds.add(n, c1)
					atoms.add(c2)
					bonds.add(n, c2)
				}

				val protonationN = mol.findProtonation(n) { it.numH == 2 && it.hybridization == Hybridization.Sp3 }
				mol.protonate(n, protonationN)
				mol.numH(n) shouldBe protonationN.numH

				val protonationC = mol.findProtonation(c1) { it.numH == 3 }
				mol.protonate(c1, protonationC)
				mol.numH(c1) shouldBe protonationC.numH
				mol.protonate(c2, protonationC)
				mol.numH(c2) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for tetrahedral geometry?
			}
		}

		test("trimethylammonium") {
			withService {

				val n = Atom(Element.Nitrogen, "N", -0.1129, -0.1597, -0.2766)
				val c1 = Atom(Element.Carbon, "C1", 1.3719, -0.2072, -0.3588)
				val c2 = Atom(Element.Carbon, "C2", -0.6526, 1.2244, -0.3588)
				val c3 = Atom(Element.Carbon, "C3", -0.6526, -0.9229, 0.8809)
				val mol = Molecule("trimethylammonium").apply {
					atoms.add(n)
					atoms.add(c1)
					bonds.add(n, c1)
					atoms.add(c2)
					bonds.add(n, c2)
					atoms.add(c3)
					bonds.add(n, c3)
				}

				val protonationN = mol.findProtonation(n) { it.numH == 1 && it.hybridization == Hybridization.Sp3 }
				mol.protonate(n, protonationN)
				mol.numH(n) shouldBe protonationN.numH

				val protonationC = mol.findProtonation(c1) { it.numH == 3 }
				mol.protonate(c1, protonationC)
				mol.numH(c1) shouldBe protonationC.numH
				mol.protonate(c2, protonationC)
				mol.numH(c2) shouldBe protonationC.numH
				mol.protonate(c3, protonationC)
				mol.numH(c3) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for tetrahedral geometry?
			}
		}
	}

	context("oxygen") {

		test("hydroxide anion") {
			withService {

				val o = Atom(Element.Oxygen, "O", 0.0, 0.0, 0.0)
				val mol = Molecule("oxygen").apply {
					atoms.add(o)
				}

				val protonation = mol.findProtonation(o) { it.numH == 1 }
				mol.protonate(o, protonation)
				mol.numH(o) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()
			}
		}

		test("water") {
			withService {

				val o = Atom(Element.Oxygen, "O", 0.0, 0.0, 0.0)
				val mol = Molecule("water").apply {
					atoms.add(o)
				}

				val protonation = mol.findProtonation(o) { it.numH == 2 }
				mol.protonate(o, protonation)
				mol.numH(o) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check H-O-H bond angle?
			}
		}

		test("hydronium") {
			withService {

				val o = Atom(Element.Oxygen, "O", 0.0, 0.0, 0.0)
				val mol = Molecule("hydronium").apply {
					atoms.add(o)
				}

				val protonation = mol.findProtonation(o) { it.numH == 3 }
				mol.protonate(o, protonation)
				mol.numH(o) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check trigonal pyrimidal geometry?
			}
		}

		test("hydrogen peroxide") {
			withService {

				val o1 = Atom(Element.Oxygen, "O1", 0.7345, 0.0, 0.0)
				val o2 = Atom(Element.Oxygen, "O2", -0.7345, 0.0, 0.0)
				val mol = Molecule("hydrogen peroxide").apply {
					atoms.add(o1)
					atoms.add(o2)
					bonds.add(o1, o2)
				}

				val protonation = mol.findProtonation(o1) { it.numH == 1 }
				mol.protonate(o1, protonation)
				mol.numH(o1) shouldBe protonation.numH
				mol.protonate(o2, protonation)
				mol.numH(o2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check tetrahedral geometry?
			}
		}

		test("methyloxonium cation") {
			withService {

				val o = Atom(Element.Oxygen, "O", 0.7100, 0.0, 0.0)
				val c = Atom(Element.Carbon, "C", -0.7100, 0.0, 0.0)
				val mol = Molecule("methyloxonium cation").apply {
					atoms.add(o)
					atoms.add(c)
					bonds.add(o, c)
				}

				val protonationO = mol.findProtonation(o) { it.numH == 2 }
				mol.protonate(o, protonationO)
				mol.numH(o) shouldBe protonationO.numH

				val protonationC = mol.findProtonation(c) { it.numH == 3 }
				mol.protonate(c, protonationC)
				mol.numH(c) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check tetrahedral geometry?
			}
		}
	}

	context("phosphorus") {

		test("phosphanide anion") {
			withService {

				val p = Atom(Element.Phosphorus, "P", 0.0, 0.0, 0.0)
				val mol = Molecule("phosphorus").apply {
					atoms.add(p)
				}

				val protonation = mol.findProtonation(p) { it.numH == 2 }
				mol.protonate(p, protonation)
				mol.numH(p) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check H-P-H angle?
			}
		}

		test("phosphine") {
			withService {

				val p = Atom(Element.Phosphorus, "P", 0.0, 0.0, 0.0)
				val mol = Molecule("phosphine").apply {
					atoms.add(p)
				}

				val protonation = mol.findProtonation(p) { it.numH == 3 }
				mol.protonate(p, protonation)
				mol.numH(p) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()
			}
		}

		test("phosphonium cation") {
			withService {

				val p = Atom(Element.Phosphorus, "P", 0.0, 0.0, 0.0)
				val mol = Molecule("phosphonium cation").apply {
					atoms.add(p)
				}

				val protonation = mol.findProtonation(p) { it.numH == 4 }
				mol.protonate(p, protonation)
				mol.numH(p) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()
			}
		}

		test("methylenephosphine") {
			withService {

				val p = Atom(Element.Phosphorus, "P", 0.9275, 0.0, 0.0)
				val c = Atom(Element.Carbon, "C", -0.9275, 0.0, 0.0)
				val mol = Molecule("methylenephosphine").apply {
					atoms.add(p)
					atoms.add(c)
					bonds.add(p, c)
				}

				val protonationP = mol.findProtonation(p) { it.numH == 1 }
				mol.protonate(p, protonationP)
				mol.numH(p) shouldBe protonationP.numH

				val protonationC = mol.findProtonation(c) { it.numH == 2 }
				mol.protonate(c, protonationC)
				mol.numH(c) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for planar geometry?
			}
		}

		// TODO: can't find amber types that make PH2+ planar, so just disable this test for now
		xtest("methylenephosphonium cation") {
			withService {

				val p = Atom(Element.Phosphorus, "P", 0.7750, 0.0, 0.0)
				val c = Atom(Element.Carbon, "C", -0.7750, 0.0, 0.0)
				val mol = Molecule("methylenephosphonium cation").apply {
					atoms.add(p)
					atoms.add(c)
					bonds.add(p, c)
				}

				val protonationP = mol.findProtonation(p) { it.numH == 2 && it.hybridization == Hybridization.Sp2 }
				mol.protonate(p, protonationP)
				mol.numH(p) shouldBe protonationP.numH

				val protonationC = mol.findProtonation(c) { it.numH == 2 }
				mol.protonate(c, protonationC)
				mol.numH(c) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for trigonal planar geometry?
			}
		}

		// TODO: can't find amber types that make PH2+ planar, so just disable this test for now
		xtest("oxophosphonium cation") {
			withService {

				val p = Atom(Element.Phosphorus, "P", 0.7200, 0.0, 0.0)
				val o = Atom(Element.Oxygen, "O", -0.7200, 0.0, 0.0)
				val mol = Molecule("oxophosphonium cation").apply {
					atoms.add(p)
					atoms.add(o)
					bonds.add(p, o)
				}

				val protonationP = mol.findProtonation(p) { it.numH == 2 && it.hybridization == Hybridization.Sp2 }
				mol.protonate(p, protonationP)
				mol.numH(p) shouldBe protonationP.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for trigonal planar geometry?
			}
		}

		test("diphosphane") {
			withService {

				val p1 = Atom(Element.Phosphorus, "P1", 1.107, 0.0, 0.0)
				val p2 = Atom(Element.Phosphorus, "P2", -1.107, 0.0, 0.0)
				val mol = Molecule("diphosphane").apply {
					atoms.add(p1)
					atoms.add(p2)
					bonds.add(p1, p2)
				}

				val protonation = mol.findProtonation(p1) { it.numH == 2  && it.hybridization == Hybridization.Sp3 }
				mol.protonate(p1, protonation)
				mol.numH(p1) shouldBe protonation.numH
				mol.protonate(p2, protonation)
				mol.numH(p2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for trigonal pyramidal geometry?
			}
		}

		test("methylphosphonium") {
			withService {

				val p = Atom(Element.Phosphorus, "P", 0.8250, 0.0, 0.0)
				val c = Atom(Element.Carbon, "C", -0.8250, 0.0, 0.0)
				val mol = Molecule("methylphosphonium").apply {
					atoms.add(p)
					atoms.add(c)
					bonds.add(p, c)
				}

				val protonationP = mol.findProtonation(p) { it.numH == 3 }
				mol.protonate(p, protonationP)
				mol.numH(p) shouldBe protonationP.numH

				val protonationC = mol.findProtonation(c) { it.numH == 3 }
				mol.protonate(c, protonationC)
				mol.numH(c) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for tetrahedral geometry?
				// TODO: check for staggered conformation?
			}
		}

		test("dimethylphosphine") {
			withService {

				val p = Atom(Element.Phosphorus, "P", -0.5217, -0.7377, -0.2316)
				val c1 = Atom(Element.Carbon, "C1", 1.1331, -0.6643, -0.0976)
				val c2 = Atom(Element.Carbon, "C2", -1.0040, 0.8468, -0.0976)
				val mol = Molecule("dimethylphosphine").apply {
					atoms.add(p)
					atoms.add(c1)
					bonds.add(p, c1)
					atoms.add(c2)
					bonds.add(p, c2)
				}

				val protonationP = mol.findProtonation(p) { it.numH == 1 }
				mol.protonate(p, protonationP)
				mol.numH(p) shouldBe protonationP.numH

				val protonationC = mol.findProtonation(c1) { it.numH == 3 }
				mol.protonate(c1, protonationC)
				mol.numH(c1) shouldBe protonationC.numH
				mol.protonate(c2, protonationC)
				mol.numH(c2) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for trigonal pyramidal geometry?
			}
		}

		test("dimethylphosphonium") {
			withService {

				val p = Atom(Element.Phosphorus, "P", -0.5217, -0.7377, -0.2316)
				val c1 = Atom(Element.Carbon, "C1", 1.1331, -0.6643, -0.0976)
				val c2 = Atom(Element.Carbon, "C2", -1.0040, 0.8468, -0.0976)
				val mol = Molecule("dimethylphosphonium").apply {
					atoms.add(p)
					atoms.add(c1)
					bonds.add(p, c1)
					atoms.add(c2)
					bonds.add(p, c2)
				}

				val protonationP = mol.findProtonation(p) { it.numH == 2 }
				mol.protonate(p, protonationP)
				mol.numH(p) shouldBe protonationP.numH

				val protonationC = mol.findProtonation(c1) { it.numH == 3 }
				mol.protonate(c1, protonationC)
				mol.numH(c1) shouldBe protonationC.numH
				mol.protonate(c2, protonationC)
				mol.numH(c2) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for tetrahedral geometry?
			}
		}

		test("trimethylphosphonium") {
			withService {

				val p = Atom(Element.Phosphorus, "P", -0.0293, -0.0415, -0.0719)
				val c1 = Atom(Element.Carbon, "C1", 1.6207, -0.0415, -0.0719)
				val c2 = Atom(Element.Carbon, "C2", -0.5793, 1.5142, -0.0719)
				val c3 = Atom(Element.Carbon, "C3", -0.5793, -0.8193, 1.2754)
				val mol = Molecule("trimethylphosphonium").apply {
					atoms.add(p)
					atoms.add(c1)
					bonds.add(p, c1)
					atoms.add(c2)
					bonds.add(p, c2)
					atoms.add(c3)
					bonds.add(p, c3)
				}

				val protonationP = mol.findProtonation(p) { it.numH == 1 }
				mol.protonate(p, protonationP)
				mol.numH(p) shouldBe protonationP.numH

				val protonationC = mol.findProtonation(c1) { it.numH == 3 }
				mol.protonate(c1, protonationC)
				mol.numH(c1) shouldBe protonationC.numH
				mol.protonate(c2, protonationC)
				mol.numH(c2) shouldBe protonationC.numH
				mol.protonate(c3, protonationC)
				mol.numH(c3) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for tetrahedral geometry?
			}
		}
	}

	context("sulfur") {

		test("bisulfide anion") {
			withService {

				val s = Atom(Element.Sulfur, "S", 0.0, 0.0, 0.0)
				val mol = Molecule("sulfur").apply {
					atoms.add(s)
				}

				val protonation = mol.findProtonation(s) { it.numH == 1 }
				mol.protonate(s, protonation)
				mol.numH(s) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()
			}
		}

		test("hydrogen sulfide") {
			withService {

				val s = Atom(Element.Sulfur, "S", 0.0, 0.0, 0.0)
				val mol = Molecule("hydrogen sulfide").apply {
					atoms.add(s)
				}

				val protonation = mol.findProtonation(s) { it.numH == 2 }
				mol.protonate(s, protonation)
				mol.numH(s) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check H-S-H angle?
			}
		}

		test("sulfonium cation") {
			withService {

				val s = Atom(Element.Sulfur, "S", 0.0, 0.0, 0.0)
				val mol = Molecule("sulfonium cation").apply {
					atoms.add(s)
				}

				val protonation = mol.findProtonation(s) { it.numH == 3 }
				mol.protonate(s, protonation)
				mol.numH(s) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check trigonal pyramidal geometry?
			}
		}

		test("hydrogen disulfide") {
			withService {

				val s1 = Atom(Element.Sulfur, "S1", 1.0290, 0.0, 0.0)
				val s2 = Atom(Element.Sulfur, "S2", -1.0290, 0.0, 0.0)
				val mol = Molecule("hydrogen disulfide").apply {
					atoms.add(s1)
					atoms.add(s2)
					bonds.add(s1, s2)
				}

				val protonation = mol.findProtonation(s1) { it.numH == 1 }
				mol.protonate(s1, protonation)
				mol.numH(s1) shouldBe protonation.numH
				mol.protonate(s2, protonation)
				mol.numH(s2) shouldBe protonation.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check for staggered conf?
			}
		}

		test("methylsulfonium cation") {
			withService {

				val s = Atom(Element.Sulfur, "S", 0.92175, 0.0, 0.0)
				val c = Atom(Element.Carbon, "C", -0.92175, 0.0, 0.0)
				val mol = Molecule("methylsulfonium cation").apply {
					atoms.add(s)
					atoms.add(c)
					bonds.add(s, c)
				}

				val protonationS = mol.findProtonation(s) { it.numH == 2 }
				mol.protonate(s, protonationS)
				mol.numH(s) shouldBe protonationS.numH

				val protonationC = mol.findProtonation(c) { it.numH == 3 }
				mol.protonate(c, protonationC)
				mol.numH(c) shouldBe protonationC.numH

				mol.shouldHaveUniqueHydrogens()

				// TODO: check trigonal pyramidal geometry?
			}
		}
	}
})
