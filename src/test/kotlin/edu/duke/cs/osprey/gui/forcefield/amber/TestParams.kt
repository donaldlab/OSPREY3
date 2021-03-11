package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.SharedSpec
import edu.duke.cs.osprey.gui.forcefield.AtomIndex
import edu.duke.cs.osprey.gui.io.fromMol2
import edu.duke.cs.osprey.gui.io.fromMol2WithMetadata
import edu.duke.cs.osprey.gui.io.withService
import io.kotlintest.shouldBe


class TestParams : SharedSpec({

	fun Polymer.findNonTerminalResidue(type: String): Polymer.Residue {
		val residues = chains.first().residues
		return residues.subList(1, residues.size - 1)
			.find { it.type.toLowerCase() == type.toLowerCase() }
			?: throw NoSuchElementException("no non-terminal $type residue")
	}

	fun Collection<Atom>.assertType(types: AmberTypes, atomIndex: AtomIndex, atomName: String, atomType: String) {
		val atom = find { it.name == atomName }
			?: throw NoSuchElementException("can't find atom named $atomName in residue $this")
		val atomi = atomIndex.getOrThrow(atom)
		types.atomTypes[atomi] shouldBe atomType
	}

	fun String.skipFirstLine() = substring(indexOf('\n') + 1)

	fun AmberParams.assert(top: String, crd: String) {
		// skip the first line of the topology file, since it has a timestamp
		this.top.skipFirstLine() shouldBe top.skipFirstLine()
		this.crd shouldBe crd
	}

	group("1cc8") {

		group("protein") {

			val mol = Molecule.fromMol2(OspreyGui.getResourceAsString("1cc8.protein.h.amber.mol2")) as Polymer

			test("types") {
				withService {

					val atomIndex = AtomIndex(mol.atoms)
					val types = mol.calcTypesAmber(MoleculeType.Protein, atomIndex, ForcefieldName.ff96)

					// check the atom types against the amber94 atom types for amino acids

					mol.findNonTerminalResidue("ARG").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "CT")
						atoms.assertType(types, atomIndex, "HG2", "HC")
						atoms.assertType(types, atomIndex, "HG3", "HC")
						atoms.assertType(types, atomIndex, "CD", "CT")
						atoms.assertType(types, atomIndex, "HD2", "H1")
						atoms.assertType(types, atomIndex, "HD3", "H1")
						atoms.assertType(types, atomIndex, "NE", "N2")
						atoms.assertType(types, atomIndex, "HE", "H")
						atoms.assertType(types, atomIndex, "CZ", "CA")
						atoms.assertType(types, atomIndex, "NH1", "N2")
						atoms.assertType(types, atomIndex, "HH11", "H")
						atoms.assertType(types, atomIndex, "HH12", "H")
						atoms.assertType(types, atomIndex, "NH2", "N2")
						atoms.assertType(types, atomIndex, "HH21", "H")
						atoms.assertType(types, atomIndex, "HH22", "H")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("HIS").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "CC")
						atoms.assertType(types, atomIndex, "ND1", "NB")
						atoms.assertType(types, atomIndex, "CE1", "CR")
						atoms.assertType(types, atomIndex, "HE1", "H5")
						atoms.assertType(types, atomIndex, "NE2", "NA")
						atoms.assertType(types, atomIndex, "HE2", "H")
						atoms.assertType(types, atomIndex, "CD2", "CW")
						atoms.assertType(types, atomIndex, "HD2", "H4")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("LYS").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "CT")
						atoms.assertType(types, atomIndex, "HG2", "HC")
						atoms.assertType(types, atomIndex, "HG3", "HC")
						atoms.assertType(types, atomIndex, "CD", "CT")
						atoms.assertType(types, atomIndex, "HD2", "HC")
						atoms.assertType(types, atomIndex, "HD3", "HC")
						atoms.assertType(types, atomIndex, "CE", "CT")
						atoms.assertType(types, atomIndex, "HE2", "HP")
						atoms.assertType(types, atomIndex, "HE3", "HP")
						atoms.assertType(types, atomIndex, "NZ", "N3")
						atoms.assertType(types, atomIndex, "HZ1", "H")
						atoms.assertType(types, atomIndex, "HZ2", "H")
						atoms.assertType(types, atomIndex, "HZ3", "H")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("ASP").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "C")
						atoms.assertType(types, atomIndex, "OD1", "O2")
						atoms.assertType(types, atomIndex, "OD2", "O2")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("GLU").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "CT")
						atoms.assertType(types, atomIndex, "HG2", "HC")
						atoms.assertType(types, atomIndex, "HG3", "HC")
						atoms.assertType(types, atomIndex, "CD", "C")
						atoms.assertType(types, atomIndex, "OE1", "O2")
						atoms.assertType(types, atomIndex, "OE2", "O2")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("CYS").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "H1")
						atoms.assertType(types, atomIndex, "HB3", "H1")
						atoms.assertType(types, atomIndex, "SG", "SH")
						atoms.assertType(types, atomIndex, "HG", "HS")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("GLY").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA2", "H1")
						atoms.assertType(types, atomIndex, "HA3", "H1")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("PRO").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "CD", "CT")
						atoms.assertType(types, atomIndex, "HD2", "H1")
						atoms.assertType(types, atomIndex, "HD3", "H1")
						atoms.assertType(types, atomIndex, "CG", "CT")
						atoms.assertType(types, atomIndex, "HG2", "HC")
						atoms.assertType(types, atomIndex, "HG3", "HC")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("ALA").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB1", "HC")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("VAL").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB", "HC")
						atoms.assertType(types, atomIndex, "CG1", "CT")
						atoms.assertType(types, atomIndex, "HG11", "HC")
						atoms.assertType(types, atomIndex, "HG12", "HC")
						atoms.assertType(types, atomIndex, "HG13", "HC")
						atoms.assertType(types, atomIndex, "CG2", "CT")
						atoms.assertType(types, atomIndex, "HG21", "HC")
						atoms.assertType(types, atomIndex, "HG22", "HC")
						atoms.assertType(types, atomIndex, "HG23", "HC")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("ILE").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB", "HC")
						atoms.assertType(types, atomIndex, "CG2", "CT")
						atoms.assertType(types, atomIndex, "HG21", "HC")
						atoms.assertType(types, atomIndex, "HG22", "HC")
						atoms.assertType(types, atomIndex, "HG23", "HC")
						atoms.assertType(types, atomIndex, "CG1", "CT")
						atoms.assertType(types, atomIndex, "HG12", "HC")
						atoms.assertType(types, atomIndex, "HG13", "HC")
						atoms.assertType(types, atomIndex, "CD1", "CT")
						atoms.assertType(types, atomIndex, "HD11", "HC")
						atoms.assertType(types, atomIndex, "HD12", "HC")
						atoms.assertType(types, atomIndex, "HD13", "HC")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("LEU").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "CT")
						atoms.assertType(types, atomIndex, "HG", "HC")
						atoms.assertType(types, atomIndex, "CD1", "CT")
						atoms.assertType(types, atomIndex, "HD11", "HC")
						atoms.assertType(types, atomIndex, "HD12", "HC")
						atoms.assertType(types, atomIndex, "HD13", "HC")
						atoms.assertType(types, atomIndex, "CD2", "CT")
						atoms.assertType(types, atomIndex, "HD21", "HC")
						atoms.assertType(types, atomIndex, "HD22", "HC")
						atoms.assertType(types, atomIndex, "HD23", "HC")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("MET").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "CT")
						atoms.assertType(types, atomIndex, "HG2", "H1")
						atoms.assertType(types, atomIndex, "HG3", "H1")
						atoms.assertType(types, atomIndex, "SD", "S")
						atoms.assertType(types, atomIndex, "CE", "CT")
						atoms.assertType(types, atomIndex, "HE1", "H1")
						atoms.assertType(types, atomIndex, "HE2", "H1")
						atoms.assertType(types, atomIndex, "HE3", "H1")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("PHE").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "CA")
						atoms.assertType(types, atomIndex, "CD1", "CA")
						atoms.assertType(types, atomIndex, "HD1", "HA")
						atoms.assertType(types, atomIndex, "CE1", "CA")
						atoms.assertType(types, atomIndex, "HE1", "HA")
						atoms.assertType(types, atomIndex, "CZ", "CA")
						atoms.assertType(types, atomIndex, "HZ", "HA")
						atoms.assertType(types, atomIndex, "CE2", "CA")
						atoms.assertType(types, atomIndex, "HE2", "HA")
						atoms.assertType(types, atomIndex, "CD2", "CA")
						atoms.assertType(types, atomIndex, "HD2", "HA")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("TYR").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "CA")
						atoms.assertType(types, atomIndex, "CD1", "CA")
						atoms.assertType(types, atomIndex, "HD1", "HA")
						atoms.assertType(types, atomIndex, "CE1", "CA")
						atoms.assertType(types, atomIndex, "HE1", "HA")
						atoms.assertType(types, atomIndex, "CZ", "C")
						atoms.assertType(types, atomIndex, "OH", "OH")
						atoms.assertType(types, atomIndex, "HH", "HO")
						atoms.assertType(types, atomIndex, "CE2", "CA")
						atoms.assertType(types, atomIndex, "HE2", "HA")
						atoms.assertType(types, atomIndex, "CD2", "CA")
						atoms.assertType(types, atomIndex, "HD2", "HA")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					/* 1cc8 doesn't have TRP =(
					mol.findNonTerminalResidue("TRP").apply {
						atoms.checkType(types, atomIndex, "N", "N")
						atoms.checkType(types, atomIndex, "H", "H")
						atoms.checkType(types, atomIndex, "CA", "CT")
						atoms.checkType(types, atomIndex, "HA", "H1")
						atoms.checkType(types, atomIndex, "CB", "CT")
						atoms.checkType(types, atomIndex, "HB2", "HC")
						atoms.checkType(types, atomIndex, "HB3", "HC")
						atoms.checkType(types, atomIndex, "CG", "C*")
						atoms.checkType(types, atomIndex, "CD1", "CW")
						atoms.checkType(types, atomIndex, "HD1", "H4")
						atoms.checkType(types, atomIndex, "NE1", "NA")
						atoms.checkType(types, atomIndex, "HE1", "H")
						atoms.checkType(types, atomIndex, "CE2", "CN")
						atoms.checkType(types, atomIndex, "CZ2", "CA")
						atoms.checkType(types, atomIndex, "HZ2", "HA")
						atoms.checkType(types, atomIndex, "CH2", "CA")
						atoms.checkType(types, atomIndex, "HH2", "HA")
						atoms.checkType(types, atomIndex, "CZ3", "CA")
						atoms.checkType(types, atomIndex, "HZ3", "HA")
						atoms.checkType(types, atomIndex, "CE3", "CA")
						atoms.checkType(types, atomIndex, "HE3", "HA")
						atoms.checkType(types, atomIndex, "CD2", "CB")
						atoms.checkType(types, atomIndex, "C", "C")
						atoms.checkType(types, atomIndex, "O", "O")
					}
					*/
					// TODO: find a protein that has TRPs?

					mol.findNonTerminalResidue("SER").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "H1")
						atoms.assertType(types, atomIndex, "HB3", "H1")
						atoms.assertType(types, atomIndex, "OG", "OH")
						atoms.assertType(types, atomIndex, "HG", "HO")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("THR").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB", "H1")
						atoms.assertType(types, atomIndex, "CG2", "CT")
						atoms.assertType(types, atomIndex, "HG21", "HC")
						atoms.assertType(types, atomIndex, "HG22", "HC")
						atoms.assertType(types, atomIndex, "HG23", "HC")
						atoms.assertType(types, atomIndex, "OG1", "OH")
						atoms.assertType(types, atomIndex, "HG1", "HO")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("ASN").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "C")
						atoms.assertType(types, atomIndex, "OD1", "O")
						atoms.assertType(types, atomIndex, "ND2", "N")
						atoms.assertType(types, atomIndex, "HD21", "H")
						atoms.assertType(types, atomIndex, "HD22", "H")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}

					mol.findNonTerminalResidue("GLN").apply {
						atoms.assertType(types, atomIndex, "N", "N")
						atoms.assertType(types, atomIndex, "H", "H")
						atoms.assertType(types, atomIndex, "CA", "CT")
						atoms.assertType(types, atomIndex, "HA", "H1")
						atoms.assertType(types, atomIndex, "CB", "CT")
						atoms.assertType(types, atomIndex, "HB2", "HC")
						atoms.assertType(types, atomIndex, "HB3", "HC")
						atoms.assertType(types, atomIndex, "CG", "CT")
						atoms.assertType(types, atomIndex, "HG2", "HC")
						atoms.assertType(types, atomIndex, "HG3", "HC")
						atoms.assertType(types, atomIndex, "CD", "C")
						atoms.assertType(types, atomIndex, "OE1", "O")
						atoms.assertType(types, atomIndex, "NE2", "N")
						atoms.assertType(types, atomIndex, "HE21", "H")
						atoms.assertType(types, atomIndex, "HE22", "H")
						atoms.assertType(types, atomIndex, "C", "C")
						atoms.assertType(types, atomIndex, "O", "O")
					}
				}
			}

			test("mods") {
				withService {

					val atomIndex = AtomIndex(mol.atoms)
					val types = mol.calcTypesAmber(MoleculeType.Protein, atomIndex, ForcefieldName.ff96)
					val frcmod = mol.calcModsAmber(types, atomIndex)

					frcmod shouldBe null
				}
			}

			test("params") {
				withService {

					val atomIndex = AtomIndex(mol.atoms)
					val types = mol.calcTypesAmber(MoleculeType.Protein, atomIndex, ForcefieldName.ff96)
					val params = AmberMolParams(mol, atomIndex, types, emptyList()).calcParamsAmber()

					params.assert(
						OspreyGui.getResourceAsString("1cc8.protein.top"),
						OspreyGui.getResourceAsString("1cc8.protein.crd")
					)
				}
			}
		}

		group("benzamidine") {

			val (mol, metadata) = Molecule.fromMol2WithMetadata(OspreyGui.getResourceAsString("benzamidine.h.gaff2.mol2"))

			test("types") {
				withService {

					val atomIndex = AtomIndex(mol.atoms)
					val types = mol.calcTypesAmber(MoleculeType.SmallMolecule, atomIndex, ForcefieldName.gaff2)

					mol.atoms.assertType(types, atomIndex, "C1", "ca")
					mol.atoms.assertType(types, atomIndex, "C2", "ca")
					mol.atoms.assertType(types, atomIndex, "C3", "ca")
					mol.atoms.assertType(types, atomIndex, "C4", "ca")
					mol.atoms.assertType(types, atomIndex, "C5", "ca")
					mol.atoms.assertType(types, atomIndex, "C6", "ca")
					mol.atoms.assertType(types, atomIndex, "C", "ce")
					mol.atoms.assertType(types, atomIndex, "N1", "n2")
					mol.atoms.assertType(types, atomIndex, "N2", "n2")
					mol.atoms.assertType(types, atomIndex, "H10", "ha")
					mol.atoms.assertType(types, atomIndex, "H11", "ha")
					mol.atoms.assertType(types, atomIndex, "H12", "ha")
					mol.atoms.assertType(types, atomIndex, "H13", "ha")
					mol.atoms.assertType(types, atomIndex, "H14", "ha")
					mol.atoms.assertType(types, atomIndex, "H15", "hn")
					mol.atoms.assertType(types, atomIndex, "H16", "hn")
				}
			}

			test("mods") {
				withService {

					val atomIndex = AtomIndex(mol.atoms)
					val types = AmberTypes(ForcefieldName.gaff2, metadata, atomIndex)
					val frcmod = mol.calcModsAmber(types, atomIndex)

					frcmod shouldBe OspreyGui.getResourceAsString("benzamidine.h.frcmod")
				}
			}

			test("params") {
				withService {

					val atomIndex = AtomIndex(mol.atoms)
					val types = AmberTypes(ForcefieldName.gaff2, metadata, atomIndex)
					val frcmod = OspreyGui.getResourceAsString("benzamidine.h.frcmod")
					val params = AmberMolParams(mol, atomIndex, types, listOf(frcmod)).calcParamsAmber()

					params.assert(
						OspreyGui.getResourceAsString("benzamidine.top"),
						OspreyGui.getResourceAsString("benzamidine.crd")
					)
				}
			}
		}

		group("protein and benzamidine") {

			val (molProtein, metadataProtein) = Molecule.fromMol2WithMetadata(OspreyGui.getResourceAsString("1cc8.protein.h.amber.mol2"))
			val (molSmall, metadataSmall) = Molecule.fromMol2WithMetadata(OspreyGui.getResourceAsString("benzamidine.h.gaff2.mol2"))

			val atomIndexProtein = AtomIndex(molProtein.atoms)
			val atomIndexSmall = AtomIndex(molSmall.atoms)
			val typesProtein = AmberTypes(ForcefieldName.ff96, metadataProtein, atomIndexProtein)
			val typesSmall = AmberTypes(ForcefieldName.gaff2, metadataSmall, atomIndexSmall)

			test("types") {

				// spot check a few of the types
				molProtein as Polymer
				molProtein.findNonTerminalResidue("THR").apply {
					atoms.assertType(typesProtein, atomIndexProtein, "N", "N")
					atoms.assertType(typesProtein, atomIndexProtein, "CA", "CT")
					atoms.assertType(typesProtein, atomIndexProtein, "HG23", "HC")
				}

				molSmall.atoms.assertType(typesSmall, atomIndexSmall, "C6", "ca")
				molSmall.atoms.assertType(typesSmall, atomIndexSmall, "N1", "n2")
				molSmall.atoms.assertType(typesSmall, atomIndexSmall, "N2", "n2")
			}

			test("params") {
				withService {

					val params = listOf(
						AmberMolParams(molProtein, atomIndexProtein, typesProtein, emptyList()),
						AmberMolParams(molSmall, atomIndexSmall, typesSmall, listOf(OspreyGui.getResourceAsString("benzamidine.h.frcmod")))
					).calcParamsAmber()

					params.assert(
						OspreyGui.getResourceAsString("1cc8.protein.benzamidine.top"),
						OspreyGui.getResourceAsString("1cc8.protein.benzamidine.crd")
					)
				}
			}
		}
	}
})
