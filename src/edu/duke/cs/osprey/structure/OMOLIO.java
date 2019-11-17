package edu.duke.cs.osprey.structure;

import edu.duke.cs.osprey.tools.Streams;
import edu.duke.cs.osprey.tools.TomlParseException;
import edu.duke.cs.osprey.tools.TomlTools;
import org.tomlj.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Reads and writes OSPREY Molecules as TOML files
 */
public class OMOLIO {

	public static String write(Molecule mol) {

		StringBuilder buf = new StringBuilder();

		// write the name
		buf.append(String.format("name = %s\n", quote(mol.name)));

		// write the type, if we're not a polymer
		if (mol.residues.size() == 1) {
			buf.append(String.format("type = %s\n", quote(mol.residues.get(0).getType())));
		}

		buf.append("\n");

		// encode the atoms as a list of tables
		HashMap<Atom,Integer> indicesByAtom = new HashMap<>();
		buf.append("atoms = [\n");
		int i = 0;
		for (Residue res : mol.residues) {
			for (Atom atom : res.atoms) {
				double[] coords = atom.getCoords();
				indicesByAtom.put(atom, i);
				buf.append(String.format(
					"\t{ i=%5d, name=%7s, x=%12.6f, y=%12.6f, z=%12.6f, elem=%3s },\n",
					i,
					quote(atom.name),
					coords[0], coords[1], coords[2],
					quote(atom.elementType)
				));
				i += 1;
			}
		}
		buf.append("]\n");

		buf.append("\n");

		// encode the bonds as pairs of atom numbers
		buf.append("bonds = [\n");
		for (Residue res : mol.residues) {
			for (Atom a1 : res.atoms) {
				int i1 = indicesByAtom.get(a1);
				for (Atom a2 : a1.bonds) {
					int i2 = indicesByAtom.get(a2);
					if (i1 < i2) {
						buf.append(String.format("\t[%5d,%5d], # %6s - %-6s\n",
							i1, i2,
							a1.name, a2.name
						));
					}
				}
			}
		}
		buf.append("]");

		buf.append("\n");

		// encode the residues as a polymer
		buf.append("\n[polymer]\n");

		// for each unique chain
		mol.residues.stream()
			.map(res -> res.getChainId())
			.collect(Collectors.toSet())
			.forEach(chain -> {

				// get the residues for that chain
				buf.append(String.format("%s = [\n", quote(chain)));
				mol.residues.stream()
					.filter(res -> chain.equals(res.getChainId()))
					.forEach(res -> {

						// split into mainchain and sidechain
						List<Integer> atoms = res.atoms.stream()
							.map(atom -> indicesByAtom.get(atom))
							.collect(Collectors.toList());

						buf.append(String.format("\t{ id=%7s, type=%6s, atoms=[%s] },\n",
							quote(res.getPDBResNumber().substring(1)),
							quote(res.getType()),
							indicesToString(atoms)
						));
					});
				buf.append("]\n");
			});

		return buf.toString();
	}

	private static String indicesToString(List<Integer> indices) {
		return Streams.joinToString(indices, ", ", i -> Integer.toString(i));
	}

	private static String quote(Object val) {

		if (val == null) {
			return "\"\"";
		}

		String str = val.toString()
			.replace("\\", "\\\\")
			.replace("\"", "\\\"");
		return "\"" + str + "\"";
	}

	public static Molecule read(String toml) {

		TomlTable doc = TomlTools.parseOrThrow(toml);

		Molecule mol = new Molecule();

		// read the name
		mol.name = doc.getString("name");

		// read the atoms
		Map<Integer,Atom> atoms = new HashMap<>();
		Map<Integer,double[]> atomCoords = new HashMap<>();
		TomlArray atomsArray = TomlTools.getArrayOrThrow(doc, "atoms");
		for (int i=0; i<atomsArray.size(); i++) {
			TomlTable atomTable = TomlTools.getTableOrThrow(atomsArray, i);
			TomlPosition atomPos = atomsArray.inputPositionOf(i);

			int index = TomlTools.getIntOrThrow(atomTable, "i", atomPos);
			String name = TomlTools.getStringOrThrow(atomTable, "name", atomPos);
			double x = TomlTools.getDoubleOrThrow(atomTable, "x", atomPos);
			double y = TomlTools.getDoubleOrThrow(atomTable, "y", atomPos);
			double z = TomlTools.getDoubleOrThrow(atomTable, "z", atomPos);
			String element = TomlTools.getStringOrThrow(atomTable,"elem", atomPos);

			if (atoms.containsKey(index)) {
				throw new TomlParseException("duplicated atom index: " + index, atomPos);
			}
			atoms.put(index, new Atom(name, element));
			atomCoords.put(index, new double[] { x, y, z });
		}

		// read the polymer to get the residues
		TomlTable polymerTable = TomlTools.getTableOrThrow(doc, "polymer");
		for (String chainId : polymerTable.keySet()) {
			TomlPosition chainPos = polymerTable.inputPositionOf(chainId);
			TomlArray chainArray = TomlTools.getArrayOrThrow(polymerTable, chainId, chainPos);

			for (int i=0; i<chainArray.size(); i++) {
				TomlTable residueTable = TomlTools.getTableOrThrow(chainArray, i, chainPos);
				TomlPosition residuePos = chainArray.inputPositionOf(i);

				String id = TomlTools.getStringOrThrow(residueTable, "id", residuePos);
				String type = TomlTools.getStringOrThrow(residueTable, "type", residuePos);
				TomlArray resAtomsArray = TomlTools.getArrayOrThrow(residueTable, "atoms", residuePos);

				// collect all the atoms and coords
				ArrayList<Atom> resAtoms = new ArrayList<>();
				ArrayList<double[]> resCoords = new ArrayList<>();
				for (int a=0; a<resAtomsArray.size(); a++) {
					int index = TomlTools.getIntOrThrow(resAtomsArray, a);
					resAtoms.add(atoms.get(index));
					resCoords.add(atomCoords.get(index));
				}

				// build the residue "full name", eg "ASN A  23"
				String fullName = String.format("%3s%2s%4s",
					type.substring(0, Math.min(type.length(), 3)),
					chainId.charAt(0),
					id.substring(0, Math.min(id.length(), 4))
				);

				// make the residue
				Residue res = new Residue(resAtoms, resCoords, fullName, mol);
				mol.residues.add(res);
			}
		}

		// read the bonds
		TomlArray bondsArray = TomlTools.getArrayOrThrow(doc, "bonds");
		for (int i=0; i<bondsArray.size(); i++) {
			TomlArray bondArray = TomlTools.getArrayOrThrow(bondsArray, i);
			TomlPosition pos = bondsArray.inputPositionOf(i);

			if (bondArray.size() != 2 || !bondArray.containsLongs()) {
				throw new TomlParseException("bond needs two integers", pos);
			}
			int i1 = TomlTools.getIntOrThrow(bondArray, 0);
			int i2 = TomlTools.getIntOrThrow(bondArray, 1);

			Atom a1 = atoms.get(i1);
			Atom a2 = atoms.get(i2);

			a1.addBond(a2);
		}

		return mol;
	}
}
