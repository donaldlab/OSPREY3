package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.gui.forcefield.Forcefield
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.forcefield.amber.findTypeOrThrow
import edu.duke.cs.osprey.gui.io.fromOMOL
import edu.duke.cs.osprey.gui.io.withService
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.tools.FileTools


fun main() {

	val mols = Molecule.fromOMOL(FileTools.readFile("/home/jeff/dlab/osprey3/src/test/resources/1dg9.min.omol"))
	mols.find { it.findTypeOrThrow() == MoleculeType.Protein }!!
	val hepes = mols.find { it.findTypeOrThrow() == MoleculeType.SmallMolecule }!!

	val ffa = ForcefieldAnalyzer(mols)
	ffa.forcefields.add(Forcefield.Amber96.configure {
		vdwScale = 1.0
	})
	ffa.forcefields.add(Forcefield.EEF1.configure {
		scale = 1.0
	})
	ffa.netCharges[hepes]?.netCharge = 0

	val ffp = withService {
		ffa.parameterize()
	}

	for (ff in ffp.forcefields) {
		val energy = ffp.calcEnergy(ff.forcefield)
		println("${ff.forcefield}: $energy")
	}
}
