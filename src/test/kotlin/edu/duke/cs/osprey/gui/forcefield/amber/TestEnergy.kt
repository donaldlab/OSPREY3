package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.confspace.ParametricMolecule
import edu.duke.cs.osprey.confspace.Strand
import edu.duke.cs.osprey.energy.EnergyCalculator
import edu.duke.cs.osprey.energy.ResidueInteractions
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams as OldFFParams
import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.SharedSpec
import edu.duke.cs.osprey.gui.forcefield.AtomIndex
import edu.duke.cs.osprey.gui.forcefield.ForcefieldCalculator
import edu.duke.cs.osprey.gui.forcefield.ForcefieldParams
import edu.duke.cs.osprey.gui.io.fromOMOL
import edu.duke.cs.osprey.gui.io.toOspreyMol
import edu.duke.cs.osprey.gui.io.withService
import io.kotlintest.matchers.doubles.plusOrMinus
import io.kotlintest.shouldBe
import kotlinx.coroutines.runBlocking


class TestEnergy : SharedSpec({

	/**
	 * Calculate the molecule energy using the new parameterization system.
	 */
	fun Molecule.calcEnergyParameterized(): Double = withService {
		val mol = this
		return runBlocking {

			// parameterize the molecule
			val amberParams = Amber96Params()
			val atomIndex = AtomIndex(mol.atoms)
			val atomsParams = amberParams.parameterizeAtoms(mol, atomIndex, null)
			val atomPairsParams = amberParams.parameterizeAtomPairs(listOf(
				ForcefieldParams.MolInfo(0, mol, atomsParams, atomIndex)
			))

			// calculate the energy
			ForcefieldCalculator.calc(atomPairsParams, listOf(
				ForcefieldCalculator.MolInfo(0, mol, mol.atoms, atomIndex, atomsParams)
			))
		}
	}

	fun readMol(name: String) =
		Molecule.fromOMOL(OspreyGui.getResourceAsString("preppedMols/$name.omol"))[0]

	group("compared to templated ecalc") {

		/**
		 * The parameters from the newer code are a bit more precise than the older code,
		 * so the new energy values don't exactly match the old ones.
		 * A tolerance of 0.1 seems to cover the gap though.
		 */
		fun Double.shouldBeEnergy(expected: Double, epsilon: Double = 0.1) {
			this shouldBe expected.plusOrMinus(epsilon)
		}

		/**
		 * Calculate the molecule energy using osprey's current template-based energy calculator
		 */
		fun Molecule.calcEnergyTemplated(): Double {

			// match the molecule to osprey's templates
			val tmol = Strand.Builder(this.toOspreyMol()).build().mol

			// convert to a parametric molecule with no motions
			val pmol = ParametricMolecule(tmol)

			// use a complete residue interaction graph
			val inters = ResidueInteractions().apply {
				addComplete(tmol.residues)
			}

			// build the energy calculator with just the amber forcefield
			val ffparams = OldFFParams().apply {
				solvationForcefield = null
			}
			EnergyCalculator.Builder(ffparams).build().use { ecalc ->

				// finally, calculate the energy
				return ecalc.calcEnergy(pmol, inters).energy
			}
		}

		test("gly-gly") {
			val mol = readMol("gly-gly")
			mol.calcEnergyParameterized().shouldBeEnergy(mol.calcEnergyTemplated())
		}

		test("trp-trp") {
			val mol = readMol("trp-trp")
			mol.calcEnergyParameterized().shouldBeEnergy(mol.calcEnergyTemplated())
		}

		test("ser-met") {
			val mol = readMol("ser-met")
			mol.calcEnergyParameterized().shouldBeEnergy(mol.calcEnergyTemplated())
		}

		test("1cc8") {
			val mol = Molecule.fromOMOL(OspreyGui.getResourceAsString("1cc8.protein.omol"))[0]
			mol.calcEnergyParameterized().shouldBeEnergy(mol.calcEnergyTemplated())
		}
	}

	/*
		These tests are designed to catch regressions in the amber energy calculator.
		The expected values were captured when the code was last deemed to be working correctly.
	*/
	group("precision regressions") {

		fun Double.shouldBeEnergy(expected: Double, epsilon: Double = 1e-9) {
			this shouldBe expected.plusOrMinus(epsilon)
		}

		test("gly-gly") {
			readMol("gly-gly").calcEnergyParameterized().shouldBeEnergy(-2.908253272320646)
		}

		test("trp-trp") {
			readMol("trp-trp").calcEnergyParameterized().shouldBeEnergy(1.5169212741454106)
		}

		test("ser-met") {
			readMol("ser-met").calcEnergyParameterized().shouldBeEnergy(-3.513246964431621)
		}

		test("1cc8") {
			val mol = Molecule.fromOMOL(OspreyGui.getResourceAsString("1cc8.protein.omol"))[0]
			mol.calcEnergyParameterized().shouldBeEnergy(-489.08432295295387)
		}
	}
})
