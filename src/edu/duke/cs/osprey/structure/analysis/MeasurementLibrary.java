/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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

	public static class DeltaDihedralAngle extends Measurement {

		public final String a;
		public final String b;
		public final String c;
		public final String[] d;

		public DeltaDihedralAngle(String a, String b, String c, String d1, String d2) {
			super(String.format("DDihedral-%s-%s-%s-%s/%s", a, b, c, d1, d2), Space.S);
			this.a = a;
			this.b = b;
			this.c = c;
			this.d = new String[] { d1, d2 };
		}

		@Override
		public Double measure(Residue res) {

			// get the atoms, if possible
			Atom a = res.getAtomByName(this.a);
			Atom b = res.getAtomByName(this.b);
			Atom c = res.getAtomByName(this.c);
			Atom[] d = new Atom[] {
				res.getAtomByName(this.d[0]),
				res.getAtomByName(this.d[1])
			};
			if (a == null || b == null || c == null || d[0] == null || d[1] == null) {
				return null;
			}

			return Protractor.getDistDegrees(
				Protractor.measureDihedral(res.coords, a.indexInRes, b.indexInRes, c.indexInRes, d[0].indexInRes),
				Protractor.measureDihedral(res.coords, a.indexInRes, b.indexInRes, c.indexInRes, d[1].indexInRes)
			);
		}
	}

	/** returns the smallest dihedral angle distance to 180 among the d atoms */
	public static class DihedralAnglesMinDist extends Measurement {

		public final double refAngle;
		public final String a;
		public final String b;
		public final String c;
		public final String[] d;

		public DihedralAnglesMinDist(String name, double refAngle, String a, String b, String c, String ... d) {
			super(name, Space.S);
			this.refAngle = refAngle;
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
			if (a == null || b == null || c == null) {
				return null;
			}
			Atom[] d = new Atom[this.d.length];
			for (int i=0; i<this.d.length; i++) {
				Atom di = res.getAtomByName(this.d[i]);
				if (di == null) {
					return null;
				}
				d[i] = di;
			}

			// get closest dist to 180 degrees
			double minDist = Double.POSITIVE_INFINITY;
			for (Atom di : d) {
				double angle = Protractor.measureDihedral(res.coords, a.indexInRes, b.indexInRes, c.indexInRes, di.indexInRes);
				minDist = Math.min(minDist, Protractor.getDistDegrees(angle, refAngle));
			}

			return minDist;
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

	public boolean contains(String type) {
		return measurements.containsKey(type.toUpperCase());
	}

	public double[] measure(Residue res, String type) {
		List<Measurement> measurements = get(type);
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
