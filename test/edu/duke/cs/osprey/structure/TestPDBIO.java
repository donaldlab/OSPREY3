package edu.duke.cs.osprey.structure;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.Iterator;

import org.junit.Test;

import edu.duke.cs.osprey.structure.Residue.SecondaryStructure;
import edu.duke.cs.osprey.tools.FileTools;

public class TestPDBIO {
	
	@Test
	public void read1CC8() {
		
		Molecule mol = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");
		
		assertThat(mol.residues.size(), is(72));
		
		// spot check a few residues,atoms
		
		assertRes(mol.residues.get(0), "ALA A   2", 0, "A2");
		assertAtom(mol.residues.get(0).atoms.get(0), "N", "N", 14.699, 27.060, 24.044);
		assertAtom(mol.residues.get(0).atoms.get(9), "3HB", "H", 12.825, 26.532, 25.978);
		assertThat(mol.getAlternates(0).isEmpty(), is(true));
		
		assertRes(mol.residues.get(71), "LEU A  73", 71, "A73");
		assertAtom(mol.residues.get(71).atoms.get(0), "N", "N", 7.624, 25.000, 9.774);
		assertAtom(mol.residues.get(71).atoms.get(19), "OXT", "O", 5.315, 27.215, 11.392);
		assertThat(mol.getAlternates(71).isEmpty(), is(true));
		
		// check secondary structure
		assertThat(mol.residues.get(0).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(2).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(3).secondaryStruct, is(SecondaryStructure.SHEET));
		assertThat(mol.residues.get(9).secondaryStruct, is(SecondaryStructure.SHEET));
		assertThat(mol.residues.get(10).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(30).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(31).secondaryStruct, is(SecondaryStructure.SHEET));
		assertThat(mol.residues.get(37).secondaryStruct, is(SecondaryStructure.SHEET));
		assertThat(mol.residues.get(38).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(41).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(42).secondaryStruct, is(SecondaryStructure.SHEET));
		assertThat(mol.residues.get(47).secondaryStruct, is(SecondaryStructure.SHEET));
		assertThat(mol.residues.get(48).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(50).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(51).secondaryStruct, is(SecondaryStructure.HELIX));
		assertThat(mol.residues.get(61).secondaryStruct, is(SecondaryStructure.HELIX));
		assertThat(mol.residues.get(62).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(64).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(65).secondaryStruct, is(SecondaryStructure.SHEET));
		assertThat(mol.residues.get(71).secondaryStruct, is(SecondaryStructure.SHEET));
	}
	
	@Test
	public void read4NPD() {
		
		Molecule mol = PDBIO.readFile("examples/4NPD/4NPD.pdb");
		
		assertThat(mol.residues.size(), is(195));
		
		// spot check a few residues,atoms
		
		assertRes(mol.residues.get(0), "ALA A   1", 0, "A1");
		assertAtom(mol.residues.get(0).atoms.get(0), "N", "N", 0.666, 9.647, -8.772, 23.86);
		assertAtom(mol.residues.get(0).atoms.get(11), "HB3", "H", -0.025, 11.526, -7.056, 23.49);
		
		assertThat(mol.getAlternates(0).size(), is(1));
		assertRes(mol.getAlternates(0).get(0), "ALA A   1", 0, "A1");
		assertAtom(mol.getAlternates(0).get(0).atoms.get(0), "N", "N", 0.419, 9.252, -8.647, 24.12);
		assertAtom(mol.getAlternates(0).get(0).atoms.get(11), "HB3", "H", 0.169, 11.490, -7.227, 23.49);
		
		assertRes(mol.residues.get(14), "GLU A  15", 14, "A15");
		assertThat(mol.residues.get(14).atoms.size(), is(15));
		assertAtom(mol.residues.get(14).atoms.get(0), "N", "N", 12.735, 16.976, 11.694, 4.03);
		assertAtom(mol.residues.get(14).atoms.get(14), "HG3", "H", 14.983, 19.618, 13.536, 6.15);
		
		assertThat(mol.getAlternates(14).size(), is(2));
		assertRes(mol.getAlternates(14).get(0), "GLU A  15", 14, "A15");
		assertThat(mol.getAlternates(14).get(0).atoms.size(), is(15));
		assertAtom(mol.getAlternates(14).get(0).atoms.get(0), "N", "N", 12.878, 16.984, 11.750, 3.04);
		assertAtom(mol.getAlternates(14).get(0).atoms.get(14), "HG3", "H", 16.105, 19.267, 11.142, 6.15);
		assertRes(mol.getAlternates(14).get(1), "GLU A  15", 14, "A15");
		assertThat(mol.getAlternates(14).get(1).atoms.size(), is(14));
		assertAtom(mol.getAlternates(14).get(1).atoms.get(0), "N", "N", 12.889, 17.024, 11.679, 3.60);
		assertAtom(mol.getAlternates(14).get(1).atoms.get(13), "HG3", "H", 15.305, 19.571, 13.551, 6.15);
		
		assertRes(mol.residues.get(21), "LEU A  22", 21, "A22");
		assertThat(mol.residues.get(21).atoms.size(), is(19));
		assertAtom(mol.residues.get(21).atoms.get(0), "N", "N", 8.697, 23.455, 20.979, 4.67);
		assertAtom(mol.residues.get(21).atoms.get(9), "HA", "H", 7.257, 24.674, 20.449, 4.33);
		assertAtom(mol.residues.get(21).atoms.get(18), "HD23", "H", 8.354, 25.557, 18.774, 4.96);
		
		assertThat(mol.getAlternates(21).size(), is(1));
		assertThat(mol.getAlternates(21).get(0).atoms.size(), is(19));
		assertAtom(mol.getAlternates(21).get(0).atoms.get(0), "N", "N", 8.668, 23.395, 21.029, 4.77);
		assertAtom(mol.getAlternates(21).get(0).atoms.get(9), "HA", "H", 7.283, 24.675, 20.417, 4.33);
		assertAtom(mol.getAlternates(21).get(0).atoms.get(18), "HD23", "H", 8.354, 25.557, 18.774, 4.96);
		
		assertThat(mol.getAlternates(57).size(), is(1));
		assertRes(mol.getAlternates(57).get(0), "LYS A  58", 57, "A58");
		assertAtom(mol.getAlternates(57).get(0).atoms.get(0), "N", "N", 15.456, 28.995, 29.847, 8.93);
		assertAtom(mol.getAlternates(57).get(0).atoms.get(22), "HZ3", "H", 16.269, 35.709, 27.594, 18.55);

		assertRes(mol.residues.get(58), " ZN A 101", 58, "A101");
		assertAtom(mol.residues.get(58).atoms.get(0), "ZN", "Zn", 17.919, 29.932, 34.195, 5.56);
		assertThat(mol.getAlternates(58).isEmpty(), is(true));
		
		assertRes(mol.residues.get(60), "SCN A 103", 60, "A103");
		assertAtom(mol.residues.get(60).atoms.get(0), "S", "S", 15.650, 26.105, 35.922, 8.13);
		assertThat(mol.getAlternates(60).isEmpty(), is(true));
		
		assertRes(mol.residues.get(64), "HOH A 204", 64, "A204");
		assertAtom(mol.residues.get(64).atoms.get(0), "O", "O", 12.806, 22.871, -2.789, 6.18);
		assertThat(mol.getAlternates(64).isEmpty(), is(true));
		
		assertRes(mol.residues.get(164), "HOH A 304", 164, "A304");
		assertAtom(mol.residues.get(164).atoms.get(0), "O", "O", 14.926, 22.916,4.590, 7.78);
		assertThat(mol.getAlternates(164).size(), is(1));
		assertRes(mol.getAlternates(164).get(0), "HOH A 304", 164, "A304");
		assertAtom(mol.getAlternates(164).get(0).atoms.get(0), "O", "O", 15.973, 24.086, 4.831, 11.48);
		
		// check secondary structure
		assertThat(mol.residues.get(0).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(4).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(5).secondaryStruct, is(SecondaryStructure.HELIX));
		assertThat(mol.residues.get(18).secondaryStruct, is(SecondaryStructure.HELIX));
		assertThat(mol.residues.get(19).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(21).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(22).secondaryStruct, is(SecondaryStructure.HELIX));
		assertThat(mol.residues.get(36).secondaryStruct, is(SecondaryStructure.HELIX));
		assertThat(mol.residues.get(37).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(38).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(39).secondaryStruct, is(SecondaryStructure.HELIX));
		assertThat(mol.residues.get(54).secondaryStruct, is(SecondaryStructure.HELIX));
		assertThat(mol.residues.get(55).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(57).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(58).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(60).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(64).secondaryStruct, is(SecondaryStructure.LOOP));
		assertThat(mol.residues.get(164).secondaryStruct, is(SecondaryStructure.LOOP));
	}
	
	@Test
	public void readWrite1CC8() {
		assertReadWrite(FileTools.readFile("examples/1CC8/1CC8.copy.pdb"));
	}
	
	private void assertRes(Residue res, String name, int index, String resNum) {
		assertThat(res.fullName, is(name));
		assertThat(res.indexInMolecule, is(index));
		assertThat(res.getPDBResNumber(), is(resNum));
	}
	
	private void assertAtom(Atom atom, String name, String elem, double x, double y, double z) {
		assertAtom(atom, name, elem, x, y, z, 0);
	}
	
	private void assertAtom(Atom atom, String name, String elem, double x, double y, double z, double bFactor) {
		assertThat(atom.name, is(name));
		assertThat(atom.elementType, is(elem));
		double[] coords = atom.getCoords();
		assertThat(coords[0], is(x));
		assertThat(coords[1], is(y));
		assertThat(coords[2], is(z));
		assertThat(atom.BFactor, is(bFactor));
	}
	
	private void assertReadWrite(String pdbText) {
		Molecule mol = PDBIO.read(pdbText);
		String pdbText2 = PDBIO.write(mol);
		Iterator<String> lines = FileTools.parseLines(pdbText).iterator();
		Iterator<String> lines2 = FileTools.parseLines(pdbText2).iterator();
		while (lines.hasNext()) {
			assertThat(lines2.hasNext(), is(true));
			assertThat(nextImportantLine(lines2), is(nextImportantLine(lines)));
		}
	}

	private String nextImportantLine(Iterator<String> lines) {
		while (true) {
			String line = lines.next();
			if (line.startsWith("REMARK") || line.startsWith("AUTHOR")) {
				continue;
			}
			return line;
		}
	}
}
