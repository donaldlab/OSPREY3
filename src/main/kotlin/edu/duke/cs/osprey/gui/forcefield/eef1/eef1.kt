package edu.duke.cs.osprey.gui.forcefield.eef1

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule


object EEF1 {

	/**
	 * Atom types for EEF1 as defined in Table I, in the paper:
	 *
	 * Effective Energy Function for Proteins in Solution
	 * Themis Lazaridis and Martin Karplus
	 * Proteins. 1999
	 * https://doi.org/10.1002/(SICI)1097-0134(19990501)35:2%3C133::AID-PROT1%3E3.0.CO;2-N
	 *
	 * Parameters were derived under the assumption of T = 298.15 K
	 *
	 * The vdW radii are taken from the EEF1 parameter set in CHARMM 19,
	 * as defined in Appendix B, Table IV, in the paper:
	 *
	 * Simulation of activation free energies in molecular systems
	 * Eyal Neria, Stefan Fischer, and Martin Karplus
	 * J. Chem. Phys. 1996
	 * https://doi.org/10.1063/1.472061
	 *
	 * Each vdwRadius parameter is commented with the CHARMM 19 atom type
	 * from which the van der Waals radius was taken.
	 */
	enum class AtomType(
		val volume: Double,
		val dGref: Double,
		val dGfree: Double,
		/** not currently used in any calculations */
		val dHref: Double,
		/** not currently used in any calculations */
		val dCpref: Double,
		val lambda: Double,
		val vdwRadius: Double
	) {

		/** carbonyl carbon */
		C(    14.7,     0.0,    0.0,     0.0,    0.0, 3.5,    2.1 /* from C */),
		/** carbon with no hydrogens */
		CR(    8.3,   -0.89,   -1.4,    2.22,    6.9, 3.5,    2.1 /* from C */),
		/** extended aliphatic carbon with one hydrogen */
		CH1E( 23.7,  -0.187,  -0.25,   0.876,    0.0, 3.5,  2.365 /* from CH */),
		/** extended aliphatic carbon with two hydrogens */
		CH2E( 22.4,   0.372,   0.52,   -0.61,   18.6, 3.5,  2.235 /* from C(2) */),
		/** extended aliphatic carbon with three hydrogens */
		CH3E( 30.0,   1.089,    1.5,  -1.779,   35.6, 3.5,  2.165 /* from C(3) */),
		/** extended aromatic carbon with one hydrogen */
		CR1E( 18.4,   0.057,   0.08,  -0.973,    6.9, 3.5,    2.1 /* from CR */),
		/** nitrogen with one hydrogen, eg amide nitrogen */
		NH1(   4.4,   -5.95,   -8.9,  -9.059,   -8.8, 3.5,    1.6 /* from N1 */),
		/** aromatic nitrogen with no hydrogens */
		NR(    4.4,   -3.82,   -4.0,  -4.654,   -8.8, 3.5,    1.6 /* from NR */),
		/** nitrogen bound to two hydrogens */
		NH2(  11.2,   -5.45,   -7.8,  -9.028,   -7.0, 3.5,    1.6 /* from N2 */),
		/** nitrogen bound to three hydrogens */
		NH3(  11.2,   -20.0,  -20.0,   -25.0,  -18.0, 6.0,    1.6 /* from N3 */),
		/** guanidinium nitrogen */
		NC2(  11.2,   -10.0,  -10.0,   -12.0,   -7.0, 6.0,    1.6 /* from NC2 */),
		/** nitrogen with no hydrogens, eg proline nitrogen */
		N(     0.0,    -1.0,  -1.55,   -1.25,    8.8, 3.5,    1.6 /* from N */),
		/** hydroxyl oxygen */
		OH1(   10.8,  -5.92,   -6.7,  -9.264,  -11.2, 3.5,    1.6 /* from OH */),
		/** carbonyl oxygen */
		O(     10.8,  -5.33,  -5.85,  -5.787,   -8.8, 3.5,    1.6 /* from O */),
		/** carboxyl oxygen */
		OC(    10.8,  -10.0,  -10.0,   -12.0,   -9.4, 6.0,    1.6 /* from OC */),
		/** sulphur */
		S(     14.7,  -3.24,   -4.1,  -4.475,  -39.9, 3.5,   1.89 /* from S */),
		/** extended sulphur with one hydrogen */
		SH1E(  21.4,  -2.05,   -2.7,  -4.475,  -39.9, 3.5,   1.89 /* from S */);

		companion object {

			fun get(mol: Molecule, atom: Atom): AtomType? {

				val numH = mol.bonds.bondedAtoms(atom)
					.count { it.element == Element.Hydrogen }
				val numHeavy = mol.bonds.bondedAtoms(atom)
					.count { it.element != Element.Hydrogen }

				return when (atom.element) {

					Element.Carbon -> when (numH) {
						3 -> CH3E
						2 -> CH2E
						1 -> when (numHeavy) {
							3 -> CH1E
							2 -> CR1E
							else -> null
						}
						0 -> if (mol.bonds.bondedAtoms(atom).any { it.element == Element.Oxygen }) {
							C
						} else {
							CR
						}
						else -> null
					}

					Element.Nitrogen -> when (numH) {
						3 -> NH3
						2 ->
							// is this a guanidinium NH2 or a regular NH2?
							// look at the attached C, if any
							mol.bonds.bondedAtoms(atom)
								.filter { it.element == Element.Carbon }
								.takeIf { it.size == 1 }
								?.get(0)
								?.let { c ->
									when (mol.bonds.bondedAtoms(c).count { it.element == Element.Nitrogen }) {
										3 -> NC2
										else -> NH2
									}
								}
								?: NH2
						1 -> NH1
						0 -> when (numHeavy) {
							3 -> N
							2 -> NR
							else -> null
						}
						else -> null
					}

					Element.Oxygen -> when (numH) {
						1 -> OH1
						0 ->
							// is this carbonyl or carboxyl oxygen?
							// look at the attached C, if any
							mol.bonds.bondedAtoms(atom)
								.filter { it.element == Element.Carbon }
								.takeIf { it.size == 1 }
								?.get(0)
								?.let { c ->
									when (mol.bonds.bondedAtoms(c).count { it.element == Element.Oxygen }) {
										2 -> OC
										1 -> O
										else -> null
									}
								}
						else -> null
					}

					Element.Sulfur -> when (numH) {
						1 -> SH1E
						0 -> S
						else -> null
					}

					else -> null
				}
			}

			fun getOrThrow(mol: Molecule, atom: Atom): AtomType =
				get(mol, atom) ?: throw NoSuchElementException("no EEF1 atom type defined for $atom")
		}
	}
}

fun Atom.atomTypeEEF1(mol: Molecule) =
	EEF1.AtomType.get(mol, this)

fun Atom.atomTypeEEF1OrThrow(mol: Molecule) =
	EEF1.AtomType.getOrThrow(mol, this)
