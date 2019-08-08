package edu.duke.cs.osprey.structure;

import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.tools.Streams;
import edu.duke.cs.osprey.tools.UnpossibleError;
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

		buf.append(String.format("name = %s\n", quote(mol.name)));

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
						List<Integer> mainchain = res.atoms.stream()
							.filter(atom -> HardCodedResidueInfo.possibleBBAtomsLookup.contains(atom.name))
							.map(atom -> indicesByAtom.get(atom))
							.collect(Collectors.toList());
						List<Integer> sidechain = res.atoms.stream()
							.map(atom -> indicesByAtom.get(atom))
							.filter(index -> !mainchain.contains(index))
							.collect(Collectors.toList());

						buf.append(String.format("\t{ id=%7s, type=%6s, mainchain=[%s], sidechains=[[%s]] },\n",
							quote(res.getPDBResNumber().substring(1)),
							quote(res.template.name),
							indicesToString(mainchain),
							indicesToString(sidechain)
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

		TomlParseResult doc = Toml.parse(toml);
		if (doc.hasErrors()) {
			throw new ParseException("TOML parsing failure:\n"
				+ Streams.joinToString(doc.errors(), "\n", err -> err.toString())
			);
		}

		Molecule mol = new Molecule();

		// read the name
		mol.name = doc.getString("name");

		// read the atoms
		Map<Integer,Atom> atoms = new HashMap<>();
		Map<Integer,double[]> atomCoords = new HashMap<>();
		TomlArray atomsArray = doc.getArray("atoms");
		if (atomsArray == null) {
			throw new ParseException("missing atoms");
		}
		if (!atomsArray.containsTables()) {
			throw new ParseException("atoms does not contain tables", doc.inputPositionOf("atoms"));
		}
		for (int i=0; i<atomsArray.size(); i++) {
			TomlTable atomTable = atomsArray.getTable(i);
			TomlPosition pos = atomsArray.inputPositionOf(i);

			int index = getIntOrThrow(atomTable, pos, "i");
			String name = getStringOrThrow(atomTable, pos, "name");
			double x = getDoubleOrThrow(atomTable, pos, "x");
			double y = getDoubleOrThrow(atomTable, pos, "y");
			double z = getDoubleOrThrow(atomTable, pos, "z");
			String element = getStringOrThrow(atomTable, pos, "elem");

			if (atoms.containsKey(index)) {
				throw new ParseException("duplicated atom index: " + index, pos);
			}
			atoms.put(index, new Atom(name, element));
			atomCoords.put(index, new double[] { x, y, z });
		}

		// read the polymer to get the residues
		TomlTable polymerTable = doc.getTable("polymer");
		if (polymerTable == null) {
			throw new ParseException("missing polymer");
		}
		for (String chainId : polymerTable.keySet()) {
			TomlArray chainArray = polymerTable.getArray(chainId);
			if (chainArray == null) {
				throw new UnpossibleError(); // silly IDE linter doesn't know this can't happen
			}
			if (!chainArray.containsTables()) {
				throw new ParseException("chain does not contain tables", polymerTable.inputPositionOf(chainId));
			}

			for (int i=0; i<chainArray.size(); i++) {
				TomlTable residueTable = chainArray.getTable(i);
				TomlPosition pos = chainArray.inputPositionOf(i);

				String id = getStringOrThrow(residueTable, pos, "id");
				String type = getStringOrThrow(residueTable, pos, "type");
				TomlArray mainchain = getArrayOrThrow(residueTable, pos, "mainchain");
				TomlArray sidechains = getArrayOrThrow(residueTable, pos, "sidechains");

				// collect all the atoms and coords
				ArrayList<Atom> resAtoms = new ArrayList<>();
				ArrayList<double[]> resCoords = new ArrayList<>();
				if (!mainchain.containsLongs()) {
					throw new ParseException("field \"mainchain\" doesn't contain integers", pos);
				}
				for (int a=0; a<mainchain.size(); a++) {
					int index = (int)mainchain.getLong(a);
					resAtoms.add(atoms.get(index));
					resCoords.add(atomCoords.get(index));
				}
				if (!sidechains.containsArrays()) {
					throw new ParseException("field \"sidechains\" doesn't contain arrays", pos);
				}
				for (int s=0; s<sidechains.size(); s++) {
					TomlArray sidechain = sidechains.getArray(s);
					if (!sidechain.containsLongs()) {
						throw new ParseException("field \"sidechains\" doesn't contain arrays of integers", pos);
					}
					for (int a=0; a<sidechain.size(); a++) {
						int index = (int)sidechain.getLong(a);
						resAtoms.add(atoms.get(index));
						resCoords.add(atomCoords.get(index));
					}
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
		TomlArray bondsArray = doc.getArray("bonds");
		if (bondsArray == null) {
			throw new ParseException("no bonds");
		}
		if (!bondsArray.containsArrays()) {
			throw new ParseException("bonds does not contain arrays", doc.inputPositionOf("bonds"));
		}
		for (int i=0; i<bondsArray.size(); i++) {
			TomlArray bondArray = bondsArray.getArray(i);
			TomlPosition pos = bondsArray.inputPositionOf(i);

			if (bondArray.size() != 2 || !bondArray.containsLongs()) {
				throw new ParseException("bond needs two integers", pos);
			}
			int i1 = (int)bondArray.getLong(0);
			int i2 = (int)bondArray.getLong(1);

			Atom a1 = atoms.get(i1);
			Atom a2 = atoms.get(i2);

			a1.addBond(a2);
		}

		return mol;
	}

	private static String getStringOrThrow(TomlTable table, TomlPosition pos, String key) {
		String val = table.getString(key);
		if (val == null) {
			throw new ParseException("missing field \"" + key + "\", or it is not a string", pos);
		}
		return val;
	}

	private static int getIntOrThrow(TomlTable table, TomlPosition pos, String key) {
		Long val = table.getLong(key);
		if (val == null) {
			throw new ParseException("missing field \"" + key + "\", or it is not an integer", pos);
		}
		return val.intValue();
	}

	private static double getDoubleOrThrow(TomlTable table, TomlPosition pos, String key) {
		Double val = table.getDouble(key);
		if (val == null) {
			throw new ParseException("missing field \"" + key + "\", or it is not a floating-point number", pos);
		}
		return val;
	}

	private static TomlArray getArrayOrThrow(TomlTable table, TomlPosition pos, String key) {
		TomlArray val = table.getArray(key);
		if (val == null) {
			throw new ParseException("missing field \"" + key + "\", or it is not an array", pos);
		}
		return val;
	}

	public static class ParseException extends RuntimeException {

		public final TomlPosition pos;

		public ParseException(String msg, TomlPosition pos) {
			super(msg + (pos != null ? " at " + pos.toString() : ""));
			this.pos = pos;
		}

		public ParseException(String msg) {
			this(msg, null);
		}
	}
}
