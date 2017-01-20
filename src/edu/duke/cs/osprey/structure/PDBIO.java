package edu.duke.cs.osprey.structure;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.duke.cs.osprey.structure.Residue.SecondaryStructure;
import edu.duke.cs.osprey.tools.FileTools;

/**
 * this is a clean PDB reader that doesn't know anything about templates or bonds
 * its only job is to turn PDB text into a molecule representing a polymer with residues and atoms
 */
public class PDBIO {
	
	private static class ResInfo {
		
		public String name = null;
		public ArrayList<Atom> atoms = new ArrayList<>();
		public ArrayList<double[]> coords = new ArrayList<>();
		
		public void makeResidue(Molecule mol) {
			
			// regardless of the alt name, if this is the first res for this res number, make it the main coords
			// otherwise, add an alternate
			// NOTE: the alternate name name name is discarded
			
			Residue newRes = new Residue(atoms, coords, name, mol);
			Residue curRes = mol.getResByPDBResNumberOrNull(newRes.getPDBResNumber());
			if (curRes == null) {
				mol.appendResidue(newRes);
			} else {
				mol.addAlternate(mol.residues.size() - 1, newRes);
			}
		}

		public static void flush(Map<Character,ResInfo> resInfos, Molecule mol) {
			for (ResInfo resInfo : resInfos.values()) {
				resInfo.makeResidue(mol);
			}
			resInfos.clear();
		}
	}
	
	private static Pattern AtomLetter = Pattern.compile("([a-zA-Z])");
	
	public static Molecule readFile(String path) {
		return read(FileTools.readFile(path));
	}
	
	public static Molecule read(String pdbText) {
		return readAll(pdbText).get(0);
	}
	
	public static List<Molecule> readAll(String pdbText) {
		List<Molecule> mols = readMols(pdbText);
		readSecondaryStructure(mols, pdbText);
		return mols;
	}
	
	private static List<Molecule> readMols(String pdbText) {
		
		List<Molecule> mols = new ArrayList<>();
		Molecule mol = new Molecule();
		mols.add(mol);
		String curResName = null;
		
		// NOTE: tree map is important here
		// we need residues to be sorted by alt key in alphabetic order
		// A (or no alt) needs to get added before the alternates in ResInfo.flush() so the indices are correct
		Map<Character,ResInfo> resInfos = new TreeMap<>();
		
		for (String line : FileTools.parseLines(pdbText)) {
			line = padLine(line);
			
			if (isLine(line, "MODEL")) {
				
				// is this the first model?
				if (mol.residues.isEmpty()) {
					// ignore
				} else {
					// advance to the next molecule
					ResInfo.flush(resInfos, mol);
					mol = new Molecule();
					mols.add(mol);
				}
				
			} else if (isLine(line, "ATOM") || isLine(line, "HETATM")) {
				
				// eg
				//           1         2         3         4         5         6         7         8
				// 012345678901234567890123456789012345678901234567890123456789012345678901234567890
				// ATOM   1146  CB APRO A  38       5.781  17.860   0.637  0.45 12.10           C
				
				// parse the line
				int atomNum = Integer.parseInt(line.substring(6, 11).trim());
				String atomName = line.substring(12,16).trim();
				char alt = line.charAt(16);
				String resName = line.substring(17,27).trim();
				double x = Double.parseDouble(line.substring(30, 38).trim());
				double y = Double.parseDouble(line.substring(38, 46).trim());
				double z = Double.parseDouble(line.substring(46, 54).trim());
				double bFactor = Double.parseDouble(defaultVal("0", line.substring(60, 66).trim()));
				String elem = line.substring(76, 78).trim();
				
				if (elem.isEmpty()) {
					// no explicit element, infer from atom name
					Matcher matcher = AtomLetter.matcher(atomName);
					matcher.find();
					elem = matcher.group();
				}
				
				// should we start a new residue (with alts)?
				if (!resName.equals(curResName)) {
					ResInfo.flush(resInfos, mol);
					curResName = resName;
				}
				
				// get the res info
				ResInfo resInfo = resInfos.get(alt);
				if (resInfo == null) {
					resInfo = new ResInfo();
					resInfo.name = resName;
					resInfos.put(alt, resInfo);
				}
				
				// update the res info with the atom
				resInfo.atoms.add(new Atom(atomName, elem, bFactor, atomNum));
				resInfo.coords.add(new double[] { x, y, z });
			}
		}
		
		ResInfo.flush(resInfos, mol);
		return mols;
	}
	
	private static String padLine(String line) {
		
		final int Len = 80;
		if (line.length() >= Len) {
			return line;
		}
		
		StringBuilder buf = new StringBuilder(Len);
		buf.append(line);
		while (buf.length() < Len) {
			buf.append(' ');
		}
		return buf.toString();
	}

	private static boolean isLine(String line, String type) {
		return line.regionMatches(true, 0, type, 0, type.length());
	}
	
	private static String defaultVal(String defaultVal, String inVal) {
		if (inVal == null || inVal.isEmpty()) {
			return defaultVal;
		}
		return inVal;
	}
	
	private static void readSecondaryStructure(List<Molecule> mols, String pdbText) {
		
		// NOTE: by default, residues are assigned LOOP secondary structure
		
		// parse pass 2: read the helices and sheets
		for (String line : FileTools.parseLines(pdbText)) {
			line = padLine(line);
			
			if (isLine(line, "HELIX")) {
				
				// eg
				//           1         2         3         4         5         6         7         8
				// 012345678901234567890123456789012345678901234567890123456789012345678901234567890
				// HELIX    1   1 ASN A    6  LEU A   19  1                                  14    
				// HELIX    2   2 THR A   23  ASP A   37  1                                  15    
				// HELIX    3   3 VAL A   40  GLN A   55  1                                  16    

				String startResNum = line.substring(21, 25).trim();
				String stopResNum = line.substring(33, 37).trim();
				char chain = line.charAt(19);
				
				for (Molecule mol : mols) {
					for (Residue res : mol.getResRangeByPDBResNumber(startResNum, stopResNum)) {
						if (res.getChainId() == chain) {
							res.secondaryStruct = SecondaryStructure.HELIX;
						}
					}
				}
				
			} else if (isLine(line, "SHEET")) {
				
				// eg
				//           1         2         3         4         5         6         7         8
				// 012345678901234567890123456789012345678901234567890123456789012345678901234567890
				// SHEET    1   A 4 VAL A  67  LEU A  73  0
				// SHEET    2   A 4 LYS A   5  VAL A  11 -1  N  ASN A  10   O  ARG A  68
				// SHEET    3   A 4 LEU A  44  THR A  49 -1  N  THR A  49   O  LYS A   5
				// SHEET    4   A 4 VAL A  33  SER A  39 -1  N  SER A  39   O  LEU A  44

				String startResNum = line.substring(22, 26).trim();
				String stopResNum = line.substring(33, 37).trim();
				char chain = line.charAt(21);
				
				for (Molecule mol : mols) {
					for (Residue res : mol.getResRangeByPDBResNumber(startResNum, stopResNum)) {
						if (res.getChainId() == chain) {
							res.secondaryStruct = SecondaryStructure.SHEET;
						}
					}
				}
			}	
		}
	}
	
	public static String write(Molecule mol) {
		// TODO
		return null;
	}
}
