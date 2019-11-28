package edu.duke.cs.osprey.confspace.compiled.motions;

import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ContinuousMotion;
import edu.duke.cs.osprey.confspace.compiled.DegreeOfFreedom;
import edu.duke.cs.osprey.tools.Protractor;
import org.joml.Quaterniond;
import org.joml.Vector3d;

import java.util.List;


public class DihedralAngle implements ContinuousMotion {

	public static class Description implements ContinuousMotion.Description {

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
		public ContinuousMotion build(AssignedCoords coords, ConfSpace.Pos pos) {
			return new DihedralAngle(this, coords, pos.index);
		}
	}

	public final Description desc;
	public final AssignedCoords coords;
	public final int posi;

	public final double initialAngleRadians;
	public final double minAngleRadians;
	public final double maxAngleRadians;

	public DihedralAngle(Description desc, AssignedCoords coords, int posi) {

		this.desc = desc;
		this.coords = coords;
		this.posi = posi;

		// TODO: profile and optimize this

		// calculate the initial angle in radians
		Vector3d a = new Vector3d();
		Vector3d b = new Vector3d();
		Vector3d c = new Vector3d();
		Vector3d d = new Vector3d();
		getCoords(posi, desc.a, a);
		getCoords(posi, desc.b, b);
		getCoords(posi, desc.c, c);
		getCoords(posi, desc.d, d);
		this.initialAngleRadians = measureAngleRadians(a, b, c, d);

		this.minAngleRadians = Math.toRadians(desc.minDegrees);
		this.maxAngleRadians = Math.toRadians(desc.maxDegrees);
	}

	private void getCoords(int posi, int atomi, Vector3d out) {
		if (atomi >= 0) {
			// positive indices encode conformation atoms
			coords.getConfCoords(posi, atomi, out);
		} else {
			// negative indices encode static atoms
			coords.getStaticCoords(-atomi - 1, out);
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
		getCoords(posi, desc.a, a);
		getCoords(posi, desc.b, b);
		getCoords(posi, desc.c, c);
		getCoords(posi, desc.d, d);

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
		for (int i : desc.rotated) {
			coords.getConfCoords(posi, i, temp);
			temp.sub(b);
			temp.rotate(qIn);
			temp.rotate(qZ);
			temp.rotate(qOut);
			temp.add(b);
			coords.setConfCoords(posi, i, temp);
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
		public int[] modifiedPosIndices() {
			return new int[] { posi };
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
