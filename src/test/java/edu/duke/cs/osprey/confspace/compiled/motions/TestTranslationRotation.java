package edu.duke.cs.osprey.confspace.compiled.motions;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static edu.duke.cs.osprey.confspace.compiled.motions.Tools.*;

import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.CoordsList;
import edu.duke.cs.osprey.tools.FileTools;
import org.joml.Quaterniond;
import org.joml.Vector3d;
import org.junit.Test;


// rotations are hard  =/
public class TestTranslationRotation {

	private final double Epsilon = 1e-9;

	// load a conf space
	private static final ConfSpace cs = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/dipeptide.5hydrophobic.transrot.ccsx"));

	@Test
	public void x0() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);

		transrot.dofX.set(0.0);

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertThat(pos.x(), isAbsolutely(ipos.x(), Epsilon));
			assertThat(pos.y(), isAbsolutely(ipos.y(), Epsilon));
			assertThat(pos.z(), isAbsolutely(ipos.z(), Epsilon));
		});
	}

	@Test
	public void x1() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);

		transrot.dofX.set(1.0);

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
 			assertThat(pos.x(), isAbsolutely(ipos.x() + 1.0, Epsilon));
			assertThat(pos.y(), isAbsolutely(ipos.y(), Epsilon));
			assertThat(pos.z(), isAbsolutely(ipos.z(), Epsilon));
		});
	}

	@Test
	public void y1() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);

		transrot.dofY.set(1.0);

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertThat(pos.x(), isAbsolutely(ipos.x(), Epsilon));
			assertThat(pos.y(), isAbsolutely(ipos.y() + 1.0, Epsilon));
			assertThat(pos.z(), isAbsolutely(ipos.z(), Epsilon));
		});
	}

	@Test
	public void z1() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);

		transrot.dofZ.set(1.0);

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertThat(pos.x(), isAbsolutely(ipos.x(), Epsilon));
			assertThat(pos.y(), isAbsolutely(ipos.y(), Epsilon));
			assertThat(pos.z(), isAbsolutely(ipos.z() + 1.0, Epsilon));
		});
	}

	@Test
	public void x5x0() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);

		transrot.dofX.set(5.0);
		transrot.dofX.set(0.0);

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertThat(pos.x(), isAbsolutely(ipos.x(), Epsilon));
			assertThat(pos.y(), isAbsolutely(ipos.y(), Epsilon));
			assertThat(pos.z(), isAbsolutely(ipos.z(), Epsilon));
		});
	}

	@Test
	public void x5y2z4y7x1() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);

		transrot.dofX.set(5.0);
		transrot.dofY.set(2.0);
		transrot.dofZ.set(4.0);
		transrot.dofY.set(7.0);
		transrot.dofX.set(1.0);

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertThat(pos.x(), isAbsolutely(ipos.x() + 1.0, Epsilon));
			assertThat(pos.y(), isAbsolutely(ipos.y() + 7.0, Epsilon));
			assertThat(pos.z(), isAbsolutely(ipos.z() + 4.0, Epsilon));
		});
	}

	@Test
	public void psi0() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);

		transrot.dofPsi.set(Math.toRadians(0.0));

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertPos(pos, ipos, Epsilon);
		});
	}

	@Test
	public void psi5() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);
		var c = transrot.desc.centroid;

		transrot.dofPsi.set(Math.toRadians(5.0));

		var v = new Vector3d();
		var qx = new Quaterniond().rotationX(Math.toRadians(5.0));

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			v.set(ipos);
			v.sub(c);
			v.rotate(qx);
			v.add(c);
			assertPos(pos, v, Epsilon);
		});
	}

	@Test
	public void theta5() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);
		var c = transrot.desc.centroid;

		transrot.dofTheta.set(Math.toRadians(5.0));

		var v = new Vector3d();
		var qy = new Quaterniond().rotationY(Math.toRadians(5.0));

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			v.set(ipos);
			v.sub(c);
			v.rotate(qy);
			v.add(c);
			assertPos(pos, v, Epsilon);
		});
	}

	@Test
	public void phi5() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);
		var c = transrot.desc.centroid;

		transrot.dofPhi.set(Math.toRadians(5.0));

		var v = new Vector3d();
		var qz = new Quaterniond().rotationZ(Math.toRadians(5.0));

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			v.set(ipos);
			v.sub(c);
			v.rotate(qz);
			v.add(c);
			assertPos(pos, v, Epsilon);
		});
	}

	@Test
	public void psi5psi0() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);

		transrot.dofPsi.set(Math.toRadians(5.0));
		transrot.dofPsi.set(Math.toRadians(0.0));

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertPos(pos, ipos, Epsilon);
		});

	}

	@Test
	public void rotsThenIdentity() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);

		var rand = new Random(12345);
		for (int i=0; i<100; i++) {
			transrot.dofPsi.set(rand.nextDouble());
			transrot.dofTheta.set(rand.nextDouble());
			transrot.dofPhi.set(rand.nextDouble());
		}

		transrot.dofPsi.set(Math.toRadians(0.0));
		transrot.dofTheta.set(Math.toRadians(0.0));
		transrot.dofPhi.set(Math.toRadians(0.0));

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertPos(pos, ipos, Epsilon);
		});
	}

	@Test
	public void psi5theta7phi3theta4psi2() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);
		var c = transrot.desc.centroid;

		transrot.dofPsi.set(Math.toRadians(5.0));
		transrot.dofTheta.set(Math.toRadians(7.0));
		transrot.dofPhi.set(Math.toRadians(3.0));
		transrot.dofTheta.set(Math.toRadians(4.0));
		transrot.dofPsi.set(Math.toRadians(2.0));

		// should be equivalent to just a rotation by phi=2, theta=4, phi=3
		var v = new Vector3d();
		var qx = new Quaterniond().rotationX(Math.toRadians(2.0));
		var qy = new Quaterniond().rotationY(Math.toRadians(4.0));
		var qz = new Quaterniond().rotationZ(Math.toRadians(3.0));

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			v.set(ipos);
			v.sub(c);
			v.rotate(qx);
			v.rotate(qy);
			v.rotate(qz);
			v.add(c);
			assertPos(pos, v, Epsilon);
		});
	}

	@Test
	public void goNuts() {

		var assignedCoords = assignDipeptideGlyGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var transrot = findTransRots(assignedCoords).get(0);
		var c = transrot.desc.centroid;

		var rand = new Random(12345);
		for (int i=0; i<100; i++) {
			transrot.dofX.set(rand.nextDouble());
			transrot.dofY.set(rand.nextDouble());
			transrot.dofZ.set(rand.nextDouble());
			transrot.dofPsi.set(rand.nextDouble());
			transrot.dofTheta.set(rand.nextDouble());
			transrot.dofPhi.set(rand.nextDouble());
		}

		transrot.dofX.set(2.0);
		transrot.dofY.set(5.0);
		transrot.dofZ.set(9.0);
		transrot.dofPsi.set(Math.toRadians(7.0));
		transrot.dofTheta.set(Math.toRadians(3.0));
		transrot.dofPhi.set(Math.toRadians(4.0));

		var v = new Vector3d();
		var qx = new Quaterniond().rotationX(Math.toRadians(7.0));
		var qy = new Quaterniond().rotationY(Math.toRadians(3.0));
		var qz = new Quaterniond().rotationZ(Math.toRadians(4.0));
		var t = new Vector3d(2.0, 5.0, 9.0);

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			v.set(ipos);
			v.sub(c);
			v.rotate(qx);
			v.rotate(qy);
			v.rotate(qz);
			v.add(t);
			v.add(c);
			assertPos(pos, v, Epsilon);
		});
	}
	
	private static AssignedCoords assignDipeptideGlyGly() {
		
		// make a GLY-GLY conformation
		var conf = new int[] { 1, 1 };
		assertThat(cs.confType(0, conf[0]), is("GLY"));
		assertThat(cs.confType(1, conf[1]), is("GLY"));

		var assignedCoords = cs.makeCoords(conf);

		// we should have 6 Dofs for the translation/rotation
		assertThat(assignedCoords.dofs.size(), is(6));

		return assignedCoords;
	}
	
	public static List<TranslationRotation> findTransRots(AssignedCoords assignedCoords) {
		return assignedCoords.dofs.stream()
			.filter(dof -> dof instanceof TranslationRotation.Dof)
			.filter(dof -> dof == (((TranslationRotation.Dof)dof).translationRotation).dofX)
			.map(dof -> ((TranslationRotation.Dof)dof).translationRotation)
			.collect(Collectors.toList());
	}
}
