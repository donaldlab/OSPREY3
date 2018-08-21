package edu.duke.cs.osprey.structure.analysis;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

import java.util.*;


public class MeasurementLibrary {

	public static enum Space {
		R, // ie flat space
		S // ie circular space
	}

	public static abstract class Measurement {

		public final String name;
		public final Space space;

		protected Measurement(String name, Space space) {
			this.name = name;
			this.space = space;
		}

		public abstract Double measure(Residue res);
	}

	public static class DihedralAngle extends Measurement {

		public final String a;
		public final String b;
		public final String c;
		public final String d;

		public DihedralAngle(String a, String b, String c, String d) {
			this(String.format("Dihedral-%s-%s-%s-%s", a, b, c, d), a, b, c, d);
		}

		public DihedralAngle(String name, String a, String b, String c, String d) {
			super(name, Space.S);
			this.a = a;
			this.b = b;
			this.c = c;
			this.d = d;
		}

		@Override
		public Double measure(Residue res) {

			// get the atoms, if possible
			Atom a = res.getAtomByName(this.a);
			Atom b = res.getAtomByName(this.b);
			Atom c = res.getAtomByName(this.c);
			Atom d = res.getAtomByName(this.d);
			if (a == null || b == null || c == null || d == null) {
				return null;
			}

			return Protractor.measureDihedral(res.coords, a.indexInRes, b.indexInRes, c.indexInRes, d.indexInRes);
		}
	}

	/** dihedral of the ith d atom, sorted by distance to 150 degrees (eg CH3 groups) */
	public static class DihedralAngles3 extends Measurement {

		public final String a;
		public final String b;
		public final String c;
		public final String[] d;
		public final int dindex;

		public DihedralAngles3(String a, String b, String c, String d1, String d2, String d3, int dindex) {
			this(String.format("Dihedral-%s-%s-%s-%d/3", a, b, c, dindex + 1), a, b, c, d1, d2, d3, dindex);
		}

		public DihedralAngles3(String name, String a, String b, String c, String d1, String d2, String d3, int dindex) {
			super(name, Space.S);
			this.a = a;
			this.b = b;
			this.c = c;
			this.d = new String[] { d1, d2, d3 };
			this.dindex = dindex;
		}

		@Override
		public Double measure(Residue res) {

			// get the atoms, if possible
			Atom a = res.getAtomByName(this.a);
			Atom b = res.getAtomByName(this.b);
			Atom c = res.getAtomByName(this.c);
			Atom[] d = new Atom[] {
				res.getAtomByName(this.d[0]),
				res.getAtomByName(this.d[1]),
				res.getAtomByName(this.d[2])
			};
			if (a == null || b == null || c == null || d[0] == null || d[1] == null || d[2] == null) {
				return null;
			}

			List<Double> angles = new ArrayList<>();
			for (int i=0; i<3; i++) {
				angles.add(Protractor.measureDihedral(res.coords, a.indexInRes, b.indexInRes, c.indexInRes, d[i].indexInRes));
			}

			// sort by distance to 150 degrees
			angles.sort(Comparator.comparing(angle -> Protractor.getDistDegrees(angle, 150)));

			// pick the ith angle
			return angles.get(dindex);
		}
	}

	public static class TetrahedralSystem {

		public final String a;
		public final String b;
		public final String c;
		public final String d;

		public TetrahedralSystem(String a, String b, String c, String d) {
			this.a = a;
			this.b = b;
			this.c = c;
			this.d = d;
		}

		public Protractor.TetrahedralGeometry.Angles measure(Residue res) {

			// get the atoms, if possible
			Atom a = res.getAtomByName(this.a);
			Atom b = res.getAtomByName(this.b);
			Atom c = res.getAtomByName(this.c);
			Atom d = res.getAtomByName(this.d);
			if (a == null || b == null || c == null || d == null) {
				return null;
			}

			Protractor.TetrahedralGeometry t = new Protractor.TetrahedralGeometry();
			t.update(res.coords, a.indexInRes, b.indexInRes, c.indexInRes, d.indexInRes);
			return t.getAngles();
		}
	}

	public static class TetrahedralInPlaneAngle extends Measurement {

		public final TetrahedralSystem tetra;

		public TetrahedralInPlaneAngle(String a, String b, String c, String d) {
			super(String.format("TetraIP-%s", d), Space.S);
			this.tetra = new TetrahedralSystem(a, b, c, d);
		}

		@Override
		public Double measure(Residue res) {

			Protractor.TetrahedralGeometry.Angles angles = tetra.measure(res);
			if (angles == null) {
				return null;
			}

			return angles.inPlaneDegrees;
		}
	}

	public static class TetrahedralOutOfPlaneAngle extends Measurement {

		public final TetrahedralSystem tetra;

		public TetrahedralOutOfPlaneAngle(String a, String b, String c, String d) {
			super(String.format("TetraOOP-%s", d), Space.S);
			this.tetra = new TetrahedralSystem(a, b, c, d);
		}

		@Override
		public Double measure(Residue res) {

			Protractor.TetrahedralGeometry.Angles angles = tetra.measure(res);
			if (angles == null) {
				return null;
			}

			return angles.outOfPlaneDegrees;
		}
	}

	public static class TetrahedralSystems2 {

		public final TetrahedralSystem[] tetras;

		public TetrahedralSystems2(String a, String b, String c, String d1, String d2) {
			this.tetras = new TetrahedralSystem[] {
				new TetrahedralSystem(a, b, c, d1),
				new TetrahedralSystem(a, b, c, d2)
			};
		}

		public List<Protractor.TetrahedralGeometry.Angles> measure(Residue res) {

			List<Protractor.TetrahedralGeometry.Angles> angles = new ArrayList<>();
			for (int i = 0; i<2; i++) {
				Protractor.TetrahedralGeometry.Angles a = tetras[i].measure(res);
				if (a == null) {
					return null;
				}
				angles.add(a);
			}

			// sort by out-of-plane angle
			angles.sort(Comparator.comparing(p -> Protractor.getDistDegrees(p.outOfPlaneDegrees, -90)));

			return angles;
		}
	}

	/** in-plane angle of the ith d atom, sorted by out-of-plane angle distance to -90 degrees (eg CH2 groups) */
	public static class TetrahedralInPlaneAngles2 extends Measurement {

		public final TetrahedralSystems2 tetras;
		public final int dindex;

		public TetrahedralInPlaneAngles2(String a, String b, String c, String d1, String d2, int dindex) {
			super(String.format("TetraIP-%s/%s-%d/2", d1, d2, dindex + 1), Space.S);
			this.tetras = new TetrahedralSystems2(a, b, c, d1, d2);
			this.dindex = dindex;
		}

		@Override
		public Double measure(Residue res) {

			List<Protractor.TetrahedralGeometry.Angles> angles = tetras.measure(res);
			if (angles == null) {
				return null;
			}

			return angles.get(dindex).inPlaneDegrees;
		}
	}

	/** out-of-plane angle of the ith d atom, sorted by out-of-plane angle distance to -90 degrees (eg CH2 groups) */
	public static class TetrahedralOutOfPlaneAngles2 extends Measurement {

		public final TetrahedralSystems2 tetras;
		public final int dindex;

		public TetrahedralOutOfPlaneAngles2(String a, String b, String c, String d1, String d2, int dindex) {
			super(String.format("TetraOOP-%s/%s-%d/2", d1, d2, dindex), Space.S);
			this.tetras = new TetrahedralSystems2(a, b, c, d1, d2);
			this.dindex = dindex;
		}

		@Override
		public Double measure(Residue res) {

			List<Protractor.TetrahedralGeometry.Angles> angles = tetras.measure(res);
			if (angles == null) {
				return null;
			}

			return angles.get(dindex).outOfPlaneDegrees;
		}
	}


	public static class BondAngle extends Measurement {

		public final String a;
		public final String b;
		public final String c;

		public BondAngle(String a, String b, String c) {
			super(String.format("Angle-%s-%s-%s", a, b, c), Space.S);
			this.a = a;
			this.b = b;
			this.c = c;
		}

		@Override
		public Double measure(Residue res) {

			// get the atoms, if possible
			Atom a = res.getAtomByName(this.a);
			Atom b = res.getAtomByName(this.b);
			Atom c = res.getAtomByName(this.c);
			if (a == null || b == null || c == null) {
				return null;
			}

			return Protractor.measureBondAngle(res.coords, a.indexInRes, b.indexInRes, c.indexInRes);
		}
	}

	public static class BondLength extends Measurement {

		public final String a;
		public final String b;

		public BondLength(String a, String b) {
			super(String.format("Length-%s-%s", a, b), Space.R);
			this.a = a;
			this.b = b;
		}

		@Override
		public Double measure(Residue res) {

			// get the atoms, if possible
			Atom a = res.getAtomByName(this.a);
			Atom b = res.getAtomByName(this.b);
			if (a == null || b == null) {
				return null;
			}

			return Protractor.measureBondLength(res.coords, a.indexInRes, b.indexInRes);
		}
	}


	private final Map<String,List<Measurement>> measurements = new HashMap<>();

	public void add(String type, Measurement measurement) {
		get(type).add(measurement);
	}

	public void addAll(MeasurementLibrary other) {
		for (Map.Entry<String,List<Measurement>> entry : other.measurements.entrySet()) {
			String type = entry.getKey();
			for (Measurement measurement : entry.getValue()) {
				add(type, measurement);
			}
		}
	}

	public Measurement get(String type, String name) {
		return get(type).stream()
			.filter(m -> m.name.equalsIgnoreCase(name))
			.findAny()
			.orElse(null);
	}

	public List<Measurement> get(String type) {
		return measurements.computeIfAbsent(type.toUpperCase(), t -> new ArrayList<>());
	}

	public double[] measure(Residue res) {
		List<Measurement> measurements = get(res.getType());
		if (measurements.isEmpty()) {
			return null;
		}
		double[] p = new double[measurements.size()];
		for (int i=0; i<measurements.size(); i++) {
			Double value = measurements.get(i).measure(res);
			if (value == null) {
				return null;
			}
			p[i] = value;
		}
		return p;
	}
}
