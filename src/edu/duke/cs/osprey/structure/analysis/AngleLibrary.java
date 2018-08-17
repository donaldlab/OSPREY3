package edu.duke.cs.osprey.structure.analysis;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

public class AngleLibrary {

	public static abstract class Angle {

		public final String name;

		protected Angle(String name) {
			this.name = name;
		}

		public abstract Double measure(Residue res);
	}

	public static class Dihedral extends Angle {

		public final String[] atomNames;

		public Dihedral(String name, String ... atomNames) {
			super(name);
			this.atomNames = atomNames;
		}

		@Override
		public Double measure(Residue res) {

			// get the atoms, if possible
			Atom a = res.getAtomByName(atomNames[0]);
			Atom b = res.getAtomByName(atomNames[1]);
			Atom c = res.getAtomByName(atomNames[2]);
			Atom d = res.getAtomByName(atomNames[3]);
			if (a == null || b == null || c == null || d == null) {
				return null;
			}

			return Protractor.measureDihedral(res.coords, a.indexInRes, b.indexInRes, c.indexInRes, d.indexInRes);
		}
	}


	private final Map<String,Map<String,Angle>> angles = new HashMap<>();

	public void add(String type, Angle angle) {
		get(type).put(angle.name, angle);
	}

	public Angle get(String type, String name) {
		return get(type).get(name);
	}

	private Map<String,Angle> get(String type) {
		return angles.computeIfAbsent(type, t -> new LinkedHashMap<>());
	}

	public double[] measure(Residue res) {
		Collection<Angle> angles = get(res.getType()).values();
		double[] p = new double[angles.size()];
		int i=0;
		for (Angle angle : angles) {
			Double value = angle.measure(res);
			if (value == null) {
				return null;
			}
			p[i++] = value;
		}
		return p;
	}
}
