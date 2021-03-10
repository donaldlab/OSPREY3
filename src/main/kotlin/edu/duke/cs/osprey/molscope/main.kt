package edu.duke.cs.osprey.molscope

import cuchaz.kludge.tools.*
import edu.duke.cs.osprey.molscope.gui.Window
import edu.duke.cs.osprey.molscope.gui.features.slide.CameraTool
import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.molscope.view.BallAndStick
import edu.duke.cs.osprey.molscope.view.SpaceFilling


fun main() = autoCloser {

	// make an N-terminal alanine molecule
	val mol = Molecule("N-terminal Alanine").apply {

		val n = atoms.add(Atom(Element.Nitrogen, "N",   14.699, 27.060, 24.044))
		val h1 = atoms.add(Atom(Element.Hydrogen, "H1",  15.468, 27.028, 24.699))
		val h2 = atoms.add(Atom(Element.Hydrogen, "H2",  15.072, 27.114, 23.102))
		val h3 = atoms.add(Atom(Element.Hydrogen, "H3",  14.136, 27.880, 24.237))
		val ca = atoms.add(Atom(Element.Carbon,   "CA",  13.870, 25.845, 24.199))
		val ha = atoms.add(Atom(Element.Hydrogen, "HA",  14.468, 24.972, 23.937))
		val cb = atoms.add(Atom(Element.Carbon,   "CB",  13.449, 25.694, 25.672))
		val hb1 = atoms.add(Atom(Element.Hydrogen, "HB1", 12.892, 24.768, 25.807))
		val hb2 = atoms.add(Atom(Element.Hydrogen, "HB2", 14.334, 25.662, 26.307))
		val hb3 = atoms.add(Atom(Element.Hydrogen, "HB3", 12.825, 26.532, 25.978))
		val c = atoms.add(Atom(Element.Carbon,   "C",   12.685, 25.887, 23.222))
		val o = atoms.add(Atom(Element.Oxygen,   "O",   11.551, 25.649, 23.607))

		bonds.add(n, h1)
		bonds.add(n, h2)
		bonds.add(n, h3)
		bonds.add(n, ca)
		bonds.add(ca, ha)
		bonds.add(ca, cb)
		bonds.add(cb, hb1)
		bonds.add(cb, hb2)
		bonds.add(cb, hb3)
		bonds.add(ca, c)
		bonds.add(c, o)
	}

	// make a dipeptide
	val dipeptide = Polymer("GLU-ILE Dipeptide").apply {

		val an = atoms.add(Atom(Element.Nitrogen, "N",   12.926, 26.240, 21.956))
		val ah = atoms.add(Atom(Element.Hydrogen, "H",   13.852, 26.300, 21.559))
		val aca = atoms.add(Atom(Element.Carbon,   "CA",  11.887, 26.268, 20.919))
		val aha = atoms.add(Atom(Element.Hydrogen, "HA",  10.976, 26.735, 21.302))
		val acb = atoms.add(Atom(Element.Carbon,   "CB",  12.419, 27.131, 19.778))
		val a2hb = atoms.add(Atom(Element.Hydrogen, "2HB", 12.708, 28.110, 20.162))
		val a3hb = atoms.add(Atom(Element.Hydrogen, "3HB", 13.313, 26.637, 19.397))
		val acg = atoms.add(Atom(Element.Carbon,   "CG",  11.406, 27.343, 18.646))
		val a2hg = atoms.add(Atom(Element.Hydrogen, "2HG", 10.922, 26.393, 18.399))
		val a3hg = atoms.add(Atom(Element.Hydrogen, "3HG", 10.622, 28.022, 18.999))
		val acd = atoms.add(Atom(Element.Carbon,   "CD",  12.088, 27.904, 17.388))
		val aoe1 = atoms.add(Atom(Element.Oxygen,   "OE1", 13.342, 27.880, 17.332))
		val aoe2 = atoms.add(Atom(Element.Oxygen,   "OE2", 11.353, 28.290, 16.460))
		val ac = atoms.add(Atom(Element.Carbon,   "C",   11.569, 24.847, 20.413))
		val ao = atoms.add(Atom(Element.Oxygen,   "O",   12.441, 24.182, 19.845))

		val bn = atoms.add(Atom(Element.Nitrogen, "N",   10.337, 24.382, 20.639))
		val bh = atoms.add(Atom(Element.Hydrogen, "H",    9.687, 24.994, 21.110))
		val bca = atoms.add(Atom(Element.Carbon,   "CA",   9.771, 23.183, 20.000))
		val bha = atoms.add(Atom(Element.Hydrogen, "HA",  10.555, 22.429, 19.908))
		val bcb = atoms.add(Atom(Element.Carbon,   "CB",   8.610, 22.575, 20.829))
		val bhb = atoms.add(Atom(Element.Hydrogen, "HB",   7.790, 23.295, 20.855))
		val bcg2 = atoms.add(Atom(Element.Carbon,   "CG2",  8.115, 21.280, 20.152))
		val b1hg2 = atoms.add(Atom(Element.Hydrogen, "1HG2", 7.230, 20.907, 20.662))
		val b2hg2 = atoms.add(Atom(Element.Hydrogen, "2HG2", 7.834, 21.470, 19.117))
		val b3hg2 = atoms.add(Atom(Element.Hydrogen, "3HG2", 8.890, 20.512, 20.180))
		val bcg1 = atoms.add(Atom(Element.Carbon,   "CG1",  9.037, 22.275, 22.287))
		val b2hg1 = atoms.add(Atom(Element.Hydrogen, "2HG1", 9.753, 21.453, 22.299))
		val b3hg1 = atoms.add(Atom(Element.Hydrogen, "3HG1", 9.527, 23.148, 22.714))
		val bcd1 = atoms.add(Atom(Element.Carbon,   "CD1",  7.864, 21.935, 23.216))
		val b1hd1 = atoms.add(Atom(Element.Hydrogen, "1HD1", 8.234, 21.813, 24.235))
		val b2hd1 = atoms.add(Atom(Element.Hydrogen, "2HD1", 7.128, 22.742, 23.201))
		val b3hd1 = atoms.add(Atom(Element.Hydrogen, "3HD1", 7.384, 21.006, 22.910))
		val bc = atoms.add(Atom(Element.Carbon,   "C",    9.313, 23.581, 18.589))
		val bo = atoms.add(Atom(Element.Oxygen,   "O",    8.222, 24.116, 18.417))

		bonds.add(an, ah)
		bonds.add(an, aca)
		bonds.add(aca, aha)
		bonds.add(aca, acb)
		bonds.add(acb, a2hb)
		bonds.add(acb, a3hb)
		bonds.add(acb, acg)
		bonds.add(acg, a2hg)
		bonds.add(acg, a3hg)
		bonds.add(acg, acd)
		bonds.add(acd, aoe1)
		bonds.add(acd, aoe2)
		bonds.add(aca, ac)
		bonds.add(ac, ao)

		bonds.add(ac, bn)

		bonds.add(bn, bh)
		bonds.add(bn, bca)
		bonds.add(bca, bha)
		bonds.add(bca, bcb)
		bonds.add(bcb, bhb)
		bonds.add(bcb, bcg2)
		bonds.add(bcg2, b1hg2)
		bonds.add(bcg2, b2hg2)
		bonds.add(bcg2, b3hg2)
		bonds.add(bcb, bcg1)
		bonds.add(bcg1, b2hg1)
		bonds.add(bcg1, b3hg1)
		bonds.add(bcg1, bcd1)
		bonds.add(bcd1, b1hd1)
		bonds.add(bcd1, b2hd1)
		bonds.add(bcd1, b3hd1)
		bonds.add(bca, bc)
		bonds.add(bc, bo)

		chains.add(Polymer.Chain("A").apply {
			residues.add(Polymer.Residue(
				"1",
				"GLU",
				listOf(an, ah, aca, aha, acb, a2hb, a3hb, acg, a2hg, a3hg, acd, aoe1, aoe2, ac, ao)
			))
			residues.add(Polymer.Residue(
				"2",
				"ILE",
				listOf(bn, bh, bca, bha, bcb, bhb, bcg2, b1hg2, b2hg2, b3hg2, bcg1, b2hg1, b3hg1, bcd1, b1hd1, b2hd1, b3hd1, bc, bo)
			))
		})
	}

	// open a window
	val win = Window(
		width = 800,
		height = 600,
		title = "MolScope"
	).autoClose()

	// prepare a slide for the molecule
	win.slides.add(Slide(mol.name).apply {
		lock { s ->

			s.views.add(SpaceFilling(mol))
			s.camera.lookAtEverything()

			s.features.menu("View") {
				add(CameraTool())
			}
		}
	})

	// prepare a slide for the dipeptide
	win.slides.add(Slide(dipeptide.name).apply {
		lock { s ->

			s.views.add(BallAndStick(dipeptide))
			s.camera.lookAtEverything()

			s.features.menu("View") {
				add(CameraTool())
			}
		}
	})

	while (win.isOpen) {
		win.render()
	}
} // end of scope here cleans up all autoClose() resources
