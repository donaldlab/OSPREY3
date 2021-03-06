package edu.duke.cs.osprey.confspace.compiled.motions;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.compiled.*;
import org.joml.Quaterniond;
import org.joml.Vector3d;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class TranslationRotation implements ContinuousMotion {

	// TODO: can optimize by skipping translation,rotation DoFs during minimization
	//  when they don't effect the energy at all
	//  eg, when we're only minimizing a single molecule

	public static class Description implements ContinuousMotion.MolDescription {

		public final double maxDistance;
		public final double maxRotationRadians;
		public final Vector3d centroid;

		public Description(
			double maxDistance,
			double maxRotationRadians,
			Vector3d centroid
		) {
			this.maxDistance = maxDistance;
			this.maxRotationRadians = maxRotationRadians;
			this.centroid = centroid;
		}

		@Override
		public ContinuousMotion build(AssignedCoords coords, int molInfoIndex) {
			return new TranslationRotation(this, coords, molInfoIndex);
		}

		@Override
		public int maxNumDofs() {
			return 6;
		}
	}

	public final Description desc;
	public final AssignedCoords coords;

	public final List<Integer> atomIndices = new ArrayList<>();
	public final Set<Integer> modifiedPosIndices = new HashSet<>();

	public final Dof dofPsi;
	public final Dof dofTheta;
	public final Dof dofPhi;
	public final Dof dofX;
	public final Dof dofY;
	public final Dof dofZ;

	private final Quaterniond rotationInverse;
	private final Vector3d translation;

	public TranslationRotation(Description desc, AssignedCoords coords, int molInfoIndex) {

		this.desc = desc;
		this.coords = coords;

		// collect all the static indices for our molecule
		for (int atomi=0; atomi<coords.confSpace.numStaticAtoms; atomi++) {
			if (coords.confSpace.staticMolInfoIndices[atomi] == molInfoIndex) {
				atomIndices.add(coords.getStaticIndex(atomi));
				modifiedPosIndices.add(PosInter.StaticPos);
			}
		}

		// collect all the conformation indices for our molecule
		for (int posi=0; posi<coords.confSpace.numPos(); posi++) {

			// get the conf
			int confi = coords.assignments[posi];
			if (confi == Conf.Unassigned) {
				continue;
			}
			ConfSpace.Conf conf = coords.confSpace.positions[posi].confs[confi];

			for (int atomi=0; atomi<conf.numAtoms; atomi++) {
				if (conf.atomMolInfoIndices[atomi] == molInfoIndex) {
					atomIndices.add(coords.getConfIndex(posi, atomi));
					modifiedPosIndices.add(posi);
				}
			}
		}

		// init state to the identity transformation
		rotationInverse = new Quaterniond();
		rotationInverse.identity();
		translation = new Vector3d(0 , 0, 0);

		// make the dofs
		// for the rotations, we'll use x-y-z Tait-Bryan angles
		double rotationStep = 0.004363323; // 0.25 degrees
		double translationStep = 0.01; // angstroms
		String molName = coords.confSpace.molInfos[molInfoIndex].name;
		dofPsi = new Dof("Rotation @ " + molName + ", Psi (X)", desc.maxRotationRadians, rotationStep);
		dofTheta = new Dof("Rotation @ " + molName + ", Theta (Y)", desc.maxRotationRadians, rotationStep);
		dofPhi = new Dof("Rotation @ " + molName + ", Phi (Z)", desc.maxRotationRadians, rotationStep);
		dofX = new Dof("Translation @ " + molName + ", X", desc.maxDistance, translationStep);
		dofY = new Dof("Translation @ " + molName + ", Y", desc.maxDistance, translationStep);
		dofZ = new Dof("Translation @ " + molName + ", Z", desc.maxDistance, translationStep);
	}

	private void apply() {

		Vector3d pos = new Vector3d();
		Quaterniond q = new Quaterniond()
			.rotationXYZ(dofPsi.value, dofTheta.value, dofPhi.value);
		Vector3d t = new Vector3d(dofX.value, dofY.value, dofZ.value);

		// transform each atom
		for (int atomi : atomIndices) {
			coords.coords.get(atomi, pos);

			pos.sub(desc.centroid);

			// undo the previous transformation
			pos.sub(translation);
			pos.rotate(rotationInverse);

			// apply the new transformation
			pos.rotate(q);
			pos.add(t);

			pos.add(desc.centroid);

			coords.coords.set(atomi, pos);
		}

		// update state
		translation.set(t);
		rotationInverse.set(q).invert();
	}


	public class Dof implements DegreeOfFreedom {

		public final TranslationRotation translationRotation = TranslationRotation.this;

		public final String name;
		public final double min;
		public final double max;
		public final double step;

		private double value = 0.0;

		public Dof(String name, double dist, double step) {
			this.name = name;
			this.min = -dist;
			this.max = dist;
			this.step = step;
		}

		@Override
		public String name() {
			return name;
		}

		@Override
		public double min() {
			return min;
		}

		@Override
		public double max() {
			return max;
		}

		@Override
		public double get() {
			return value;
		}

		@Override
		public void set(double val) {
			value = val;
			apply();
		}

		@Override
		public Set<Integer> modifiedPosIndices() {
			return modifiedPosIndices;
		}

		@Override
		public double initialStepSize() {
			return step;
		}
	}

	@Override
	public void appendDofs(List<DegreeOfFreedom> dofs) {
		dofs.add(dofPsi);
		dofs.add(dofTheta);
		dofs.add(dofPhi);
		dofs.add(dofX);
		dofs.add(dofY);
		dofs.add(dofZ);
	}
}
