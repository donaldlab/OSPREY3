package edu.duke.cs.osprey.confspace.compiled.motions;

import edu.duke.cs.osprey.confspace.compiled.*;
import edu.duke.cs.osprey.tools.Protractor;
import org.joml.Quaterniond;
import org.joml.Vector3d;

import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class DihedralAngle implements ContinuousMotion {

	public static class Description implements ContinuousMotion.ConfDescription, ContinuousMotion.MolDescription {

		public final double minDegrees;
		public final double maxDegrees;

		public final int a;
		public final int b;
		public final int c;
		public final int d;
		public final int[] rotated;

		public Description(
			double minDegrees,
			double maxDegrees,
			int a,
			int b,
			int c,
			int d,
			int[] rotated
		) {
			this.minDegrees = minDegrees;
			this.maxDegrees = maxDegrees;
			this.a = a;
			this.b = b;
			this.c = c;
			this.d = d;
			this.rotated = rotated;
		}

		@Override
		public ContinuousMotion build(AssignedCoords conf, ConfSpace.Pos pos) {
			return new DihedralAngle(this, conf, 0, pos.index);
		}

		@Override
		public ContinuousMotion build(AssignedCoords conf, int molInfoIndex) {
			return new DihedralAngle(this, conf, molInfoIndex, PosInter.StaticPos);
		}

		@Override
		public int maxNumDofs() {
			return 1;
		}
	}

	public final Description desc;
	public final AssignedCoords coords;
	public final int moli;
	public final int posi;

	private final int ai;
	private final int bi;
	private final int ci;
	private final int di;
	private final int []ri;

	public final double initialAngleRadians;
	public final double minAngleRadians;
	public final double maxAngleRadians;

	public DihedralAngle(Description desc, AssignedCoords coords, int moli, int posi) {

		this.desc = desc;
		this.coords = coords;
		this.moli = moli;
		this.posi = posi;

		// cache all the atom indices
		ai = getAtomIndex(desc.a);
		bi = getAtomIndex(desc.b);
		ci = getAtomIndex(desc.c);
		di = getAtomIndex(desc.d);
		ri = new int[desc.rotated.length];
		for (int i=0; i<ri.length; i++) {
			ri[i] = getAtomIndex(desc.rotated[i]);
		}

		// TODO: profile and optimize this

		// calculate the initial angle in radians
		Vector3d a = new Vector3d();
		Vector3d b = new Vector3d();
		Vector3d c = new Vector3d();
		Vector3d d = new Vector3d();
		coords.coords.get(ai, a);
		coords.coords.get(bi, b);
		coords.coords.get(ci, c);
		coords.coords.get(di, d);
		this.initialAngleRadians = measureAngleRadians(a, b, c, d);

		this.minAngleRadians = Math.toRadians(desc.minDegrees);
		this.maxAngleRadians = Math.toRadians(desc.maxDegrees);
	}

	private String getName() {
		if (posi == PosInter.StaticPos) {
			return coords.confSpace.molInfos[moli].name;
		} else {
			return coords.confSpace.name(posi);
		}
	}

	private int getAtomIndex(int atomi) {
		if (posi == PosInter.StaticPos) {
			return coords.getStaticIndex(atomi);
		} else {
			if (atomi >= 0) {
				// positive indices encode conformation atoms
				return coords.getConfIndex(posi, atomi);
			} else {
				// negative indices encode static atoms
				return coords.getStaticIndex(-atomi - 1);
			}
		}
	}

	private String getAtomName(int atomi) {
		if (posi == PosInter.StaticPos) {
			return coords.confSpace.staticNames[atomi];
		} else {
			if (atomi >= 0) {
				// positive indices encode conformation atoms
				return coords.confSpace.positions[posi].confs[coords.assignments[posi]].atomNames[atomi];
			} else {
				// negative indices encode static atoms
				return coords.confSpace.staticNames[-atomi - 1];
			}
		}
	}

	public void setAngle(double angleRadians) {

		// TODO: profile and optimize this

		Vector3d temp = new Vector3d();
		Quaterniond qIn = new Quaterniond();
		Quaterniond qZ = new Quaterniond();
		Quaterniond qOut = new Quaterniond();

		Vector3d a = new Vector3d();
		Vector3d b = new Vector3d();
		Vector3d c = new Vector3d();
		Vector3d d = new Vector3d();

		// copy our a,b,c,d from the coords array
		coords.coords.get(ai, a);
		coords.coords.get(bi, b);
		coords.coords.get(ci, c);
		coords.coords.get(di, d);

		// translate so b is at the origin
		a.sub(b);
		c.sub(b);
		d.sub(b);

		// rotate into a coordinate system where:
		//   b->c is along the -z axis
		//   b->a is in the yz plane
		qIn.lookAlong(c, a);
		d.rotate(qIn);

		// rotate about z to set the desired dihedral angle
		qZ.rotationZ(Math.PI/2 - angleRadians - Math.atan2(d.y, d.x));

		// rotate back into the world frame
		qOut.set(qIn)
			.conjugate();

		// transform all the rotated atoms
		for (int i : ri) {
			coords.coords.get(i, temp);
			temp.sub(b);
			temp.rotate(qIn);
			temp.rotate(qZ);
			temp.rotate(qOut);
			temp.add(b);
			coords.coords.set(i, temp);
		}
	}

	/**
	 * Measures the dihedral angle in radians, using the usual a,b,c,d atom position convention.
	 *
	 * WARNING: this implementation overwrites the values of a, c, and d,
	 * but it doesn't do any heap allocations.
	 */
	private static double measureAngleRadians(Vector3d a, Vector3d b, Vector3d c, Vector3d d) {

		// translate so b is at the origin
		a.sub(b);
		c.sub(b);
		d.sub(b);

		// rotate into a coordinate system where:
		//   b->c is along the -z axis
		//   b->a is in the yz plane
		Quaterniond q = new Quaterniond()
			.lookAlong(c, a);
		d.rotate(q);

		return Protractor.normalizeMinusPiToPi(Math.PI/2 - Math.atan2(d.y, d.x));
	}

	/** a dihedral angle only has one degree of freedom */
	public class Dof implements DegreeOfFreedom {

		private double angleRadians = initialAngleRadians;

		private final Set<Integer> modifiedPosIndices = new HashSet<>();
		{
			modifiedPosIndices.add(posi);
		}

		@Override
		public String name() {
			return String.format("dihedral angle @ %s: %s-%s-%s-%s",
				getName(),
				getAtomName(desc.a),
				getAtomName(desc.b),
				getAtomName(desc.c),
				getAtomName(desc.d)
			);
		}

		@Override
		public double min() {
			return minAngleRadians;
		}

		@Override
		public double max() {
			return maxAngleRadians;
		}

		@Override
		public double get() {
			return angleRadians;
		}

		@Override
		public void set(double val) {
			if (val != angleRadians) {
				DihedralAngle.this.setAngle(val);
				angleRadians = val;
			}
		}

		@Override
		public Set<Integer> modifiedPosIndices() {
			return modifiedPosIndices;
		}

		@Override
		public double initialStepSize() {
			return 0.004363323; // 0.25 degrees
		}
	}

	@Override
	public void appendDofs(List<DegreeOfFreedom> dofs) {
		dofs.add(new Dof());
	}
}
