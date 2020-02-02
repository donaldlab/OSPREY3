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

	// TODO: is this actually correct when these Dofs get applied with other Dofs like dihedrals???
	//   maybe?? test it!!

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
	}

	public final Description desc;
	public final AssignedCoords coords;

	private final List<Integer> atomIndices = new ArrayList<>();
	private final Set<Integer> modifiedPosIndices = new HashSet<>();

	private final CoordsList originalCoords;

	private final Dof dofPsi;
	private final Dof dofTheta;
	private final Dof dofPhi;
	private final Dof dofX;
	private final Dof dofY;
	private final Dof dofZ;

	public TranslationRotation(Description desc, AssignedCoords coords, int molInfoIndex) {

		this.desc = desc;
		this.coords = coords;

		// this motion modifies the whole molecule,
		// so assume that also includes the static atoms too
		modifiedPosIndices.add(PosInter.StaticPos);

		// collect all the static indices for our molecule
		for (int atomi=0; atomi<coords.confSpace.numStaticAtoms; atomi++) {
			if (coords.confSpace.staticMolInfoIndices[atomi] == molInfoIndex) {
				atomIndices.add(coords.getStaticIndex(atomi));
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

		// copy all the coordinates for our absolute referencing
		originalCoords = new CoordsList(coords.coords.size);
		originalCoords.copyFrom(coords.coords, 0);

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
		Quaterniond qPsi = new Quaterniond().rotationX(dofPsi.value);
		Quaterniond qTheta = new Quaterniond().rotationY(dofTheta.value);
		Quaterniond qPhi = new Quaterniond().rotationZ(dofPhi.value);
		Vector3d t = new Vector3d(dofX.value, dofY.value, dofZ.value);

		// transform each atom
		for (int atomi : atomIndices) {
			originalCoords.get(atomi, pos);
			pos.sub(desc.centroid);
			pos.rotate(qPsi);
			pos.rotate(qTheta);
			pos.rotate(qPhi);
			pos.add(desc.centroid);
			pos.add(t);
			coords.coords.set(atomi, pos);
		}
	}


	public class Dof implements DegreeOfFreedom {

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
	public boolean isAbsolute() {
		// these transformations assume the atom coords start in their original positions,
		// as created by the conformation space, so we're an "absolute" transformation
		return true;
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
