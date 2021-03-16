package edu.duke.cs.osprey.confspace.compiled.motions;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import static edu.duke.cs.osprey.confspace.compiled.motions.Tools.*;
import static edu.duke.cs.osprey.confspace.compiled.motions.TestDihedralAngle.*;
import static edu.duke.cs.osprey.confspace.compiled.motions.TestTranslationRotation.*;

import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.CoordsList;
import edu.duke.cs.osprey.tools.FileTools;
import org.joml.Quaterniond;
import org.joml.Vector3d;
import org.junit.Test;

import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestCombos {

	private final double Epsilon = 1e-9;

	// load a conf space
	private static final ConfSpace cs = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/dipeptide.5hydrophobic.transrot.ccsx"));

	@Test
	public void init() {

		var assignedCoords = assignDipeptideValGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		assertThat(findDihedrals(assignedCoords).size(), is(1));
		assertThat(findTransRots(assignedCoords).size(), is(1));

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertPos(pos, ipos, Epsilon);
		});
	}

	@Test
	public void nop() {

		var assignedCoords = assignDipeptideValGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var dihedral = findDihedrals(assignedCoords).get(0);
		var transrot = findTransRots(assignedCoords).get(0);

		dihedral.setAngle(Math.toRadians(measureDihedralDegrees(dihedral)));
		transrot.dofX.set(0.0);
		transrot.dofY.set(0.0);
		transrot.dofZ.set(0.0);
		transrot.dofPsi.set(0.0);
		transrot.dofTheta.set(0.0);
		transrot.dofPhi.set(0.0);

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			assertPos(pos, ipos, Epsilon);
		});
	}

	@Test
	public void onlyDihedral() {

		var assignedCoords = assignDipeptideValGly();
		var dihedral = findDihedrals(assignedCoords).get(0);
		var transrot = findTransRots(assignedCoords).get(0);

		dihedral.setAngle(Math.toRadians(5.0));
		transrot.dofX.set(0.0);
		transrot.dofY.set(0.0);
		transrot.dofZ.set(0.0);
		transrot.dofPsi.set(0.0);
		transrot.dofTheta.set(0.0);
		transrot.dofPhi.set(0.0);

		assertThat(measureDihedralDegrees(dihedral), isAbsolutely(5.0, Epsilon));
	}

	@Test
	public void onlyTransRot() {

		var assignedCoords = assignDipeptideValGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var dihedral = findDihedrals(assignedCoords).get(0);
		var transrot = findTransRots(assignedCoords).get(0);
		var c = transrot.desc.centroid;

		dihedral.setAngle(Math.toRadians(measureDihedralDegrees(dihedral)));
		transrot.dofX.set(1.0);
		transrot.dofY.set(2.0);
		transrot.dofZ.set(3.0);
		transrot.dofPsi.set(Math.toRadians(4.0));
		transrot.dofTheta.set(Math.toRadians(5.0));
		transrot.dofPhi.set(Math.toRadians(6.0));

		var v = new Vector3d();
		var q = new Quaterniond().rotationXYZ(Math.toRadians(4.0), Math.toRadians(5.0), Math.toRadians(6.0));
		var t = new Vector3d(1.0, 2.0, 3.0);

		forEachAtom(assignedCoords, initialCoords, (pos, ipos) -> {
			v.set(ipos);
			v.sub(c);
			v.rotate(q);
			v.add(t);
			v.add(c);
			assertPos(pos, v, Epsilon);
		});
	}

	@Test
	public void dihedralAndTransRot() {

		var assignedCoords = assignDipeptideValGly();
		var initialCoords = new CoordsList(assignedCoords.coords);
		var dihedral = findDihedrals(assignedCoords).get(0);
		var transrot = findTransRots(assignedCoords).get(0);
		var c = transrot.desc.centroid;
		var initialDegrees = dihedral.measureAngleDegrees();

		transrot.dofX.set(16.0);
		transrot.dofY.set(27.0);
		transrot.dofZ.set(32.0);
		transrot.dofPsi.set(Math.toRadians(43.0));
		transrot.dofTheta.set(Math.toRadians(57.0));
		transrot.dofPhi.set(Math.toRadians(62.0));
		var deltaDegrees = 76.0;
		dihedral.setAngle(Math.toRadians(initialDegrees + deltaDegrees));

		var v = new Vector3d();
		var q = new Quaterniond().rotationXYZ(Math.toRadians(43.0), Math.toRadians(57.0), Math.toRadians(62.0));
		var t = new Vector3d(16.0, 27.0, 32.0);

		// check trans rot only atoms
		var onlyTransRotIndices = onlyTransRotIndices(assignedCoords, transrot);
		assertThat(onlyTransRotIndices.size(), is(17));
		forEachAtom(onlyTransRotIndices, assignedCoords, initialCoords, (pos, ipos) -> {
			v.set(ipos);
			v.sub(c);
			v.rotate(q);
			v.add(t);
			v.add(c);
			assertPos(pos, v, Epsilon);
		});

		// those atoms should include the dihedral a, b, and c atoms
		assertThat(onlyTransRotIndices, hasItem(dihedral.ai));
		assertThat(onlyTransRotIndices, hasItem(dihedral.bi));
		assertThat(onlyTransRotIndices, hasItem(dihedral.ci));

		// for each dihedral-rotated atom ...
		for (int ri : dihedral.ri) {

			// check the dihedral angles are offset correctly
			var obs = measureDihedralDegrees(assignedCoords.coords, dihedral, ri);
			var exp = measureDihedralDegrees(initialCoords, dihedral, ri) + deltaDegrees;
			assertThat(obs, isAbsolutely(exp, Epsilon));

			// check the r-B-C angles are preserved
			obs = measureDegrees(assignedCoords.coords, ri, dihedral.bi, dihedral.ci);
			exp = measureDegrees(initialCoords, ri, dihedral.bi, dihedral.ci);
			assertThat(obs, isAbsolutely(exp, Epsilon));

			// check the r-B distances are preserved
			obs = measureDist(assignedCoords.coords, ri, dihedral.bi);
			exp = measureDist(initialCoords, ri, dihedral.bi);
			assertThat(obs, isAbsolutely(exp, Epsilon));
		}

		// all of these three constraints should intersect to define a unique position for each atom, right?
	}

	private double measureDegrees(CoordsList coords, int ai, int bi, int ci) {

		var a = new Vector3d();
		var b = new Vector3d();
		var c = new Vector3d();

		coords.get(ai, a);
		coords.get(bi, b);
		coords.get(ci, c);

		var ba = new Vector3d(a).sub(b);
		var bc = new Vector3d(c).sub(b);

		return Math.toDegrees(Math.acos(ba.dot(bc)/ba.length()/bc.length()));
	}

	private double measureDist(CoordsList coords, int ai, int bi) {

		var a = new Vector3d();
		var b = new Vector3d();

		coords.get(ai, a);
		coords.get(bi, b);

		return a.distance(b);
	}


	private static AssignedCoords assignDipeptideValGly() {

		// make a VAL-GLY conformation
		var conf = new int[] { 14, 1 };
		assertThat(cs.confType(0, conf[0]), is("VAL"));
		assertThat(cs.confType(1, conf[1]), is("GLY"));

		var assignedCoords = cs.makeCoords(conf);

		// VAL should have just one dihedral angle, plus the 6 dofs for translations and rotations
		assertThat(assignedCoords.dofs.size(), is(7));

		return assignedCoords;
	}

	private static Set<Integer> onlyTransRotIndices(AssignedCoords assignedCoords, TranslationRotation transrot) {
		return assignedCoords.atoms().stream()
			// filter out any atoms in diherals
			.filter(i ->
				findDihedrals(assignedCoords).stream()
					.noneMatch(dihedral -> IntStream.of(dihedral.ri).anyMatch(ri -> ri == i.coordsi))
			)
			.map(i -> i.coordsi)
			.collect(Collectors.toSet());
	}
}
