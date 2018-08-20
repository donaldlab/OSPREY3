package edu.duke.cs.osprey.structure.analysis;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

import java.util.*;

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

		public Dihedral(String name, String a, String b, String c, String d) {
			super(name);
			this.atomNames = new String[] { a, b, c, d };
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

	/** eg benzene ring rotations  */
	public static class DihedralMod2 extends Dihedral {

		public DihedralMod2(String name, String a, String b, String c, String d) {
			super(name, a, b, c, d);
		}

		@Override
		public Double measure(Residue res) {

			Double angle = super.measure(res);
			if (angle == null) {
				return null;
			}

			// apply mod 2 to map to [-90,90]
			while (angle < -90) {
				angle += 180;
			}
			while (angle > 90) {
				angle -= 180;
			}

			return angle;
		}
	}

	/** eg methyl rotations  */
	public static class DihedralMod3 extends Dihedral {

		public DihedralMod3(String name, String a, String b, String c, String d) {
			super(name, a, b, c, d);
		}

		@Override
		public Double measure(Residue res) {

			Double angle = super.measure(res);
			if (angle == null) {
				return null;
			}

			// apply mod 3 to map to [-60,60]
			while (angle < -60) {
				angle += 120;
			}
			while (angle > 60) {
				angle -= 120;
			}

			return angle;
		}
	}

	public static class BondAngle extends Angle {

		public final String[] atomNames;

		public BondAngle(String a, String b, String c) {
			super(String.format("%s-%s-%s", a, b, c));
			this.atomNames = new String[] { a, b, c };
		}

		@Override
		public Double measure(Residue res) {

			// get the atoms, if possible
			Atom a = res.getAtomByName(atomNames[0]);
			Atom b = res.getAtomByName(atomNames[1]);
			Atom c = res.getAtomByName(atomNames[2]);
			if (a == null || b == null || c == null) {
				return null;
			}

			return Protractor.measureBondAngle(res.coords, a.indexInRes, b.indexInRes, c.indexInRes);
		}
	}

	public static class TetrahedralInPlaneAngle extends Angle {

		public final String[] atomNames;

		public TetrahedralInPlaneAngle(String a, String b, String c, String d) {
			super(String.format("TetraIP-%s", d));
			this.atomNames = new String[] { a, b, c, d };
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

			Protractor.TetrahedralGeometry t = new Protractor.TetrahedralGeometry();
			t.update(res.coords, a.indexInRes, b.indexInRes, c.indexInRes, d.indexInRes);
			return t.inPlaneDegrees;
		}
	}

	public static class TetrahedralOutOfPlaneAngle extends Angle {

		public final String[] atomNames;

		public TetrahedralOutOfPlaneAngle(String a, String b, String c, String d) {
			super(String.format("TetraOOP-%s", d));
			this.atomNames = new String[] { a, b, c, d };
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

			Protractor.TetrahedralGeometry t = new Protractor.TetrahedralGeometry();
			t.update(res.coords, a.indexInRes, b.indexInRes, c.indexInRes, d.indexInRes);
			return t.outOfPlaneDegrees;
		}
	}


	private final Map<String,Map<String,Angle>> angles = new HashMap<>();

	public void add(String type, Angle angle) {
		get(type).put(angle.name, angle);
	}

	public Angle get(String type, String name) {
		return get(type).get(name);
	}

	public Map<String,Angle> get(String type) {
		return angles.computeIfAbsent(type.toUpperCase(), t -> new LinkedHashMap<>());
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
