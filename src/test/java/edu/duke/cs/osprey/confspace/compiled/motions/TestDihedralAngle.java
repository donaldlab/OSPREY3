package edu.duke.cs.osprey.confspace.compiled.motions;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.isAbsolutely;

import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.CoordsList;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Protractor;
import org.joml.Vector3d;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;


public class TestDihedralAngle {

	// load a conf space
	private static final ConfSpace cs = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/dipeptide.5hydrophobic.ccsx"));

	@Test
	public void test() {

		// make a VAL-GLY conformation
		var conf = new int[] { 14, 1 };
		assertThat(cs.confType(0, conf[0]), is("VAL"));
		assertThat(cs.confType(1, conf[1]), is("GLY"));
		var coords = cs.makeCoords(conf);

		// VAL should have just one dihedral angle
		assertThat(coords.dofs.size(), is(1));
		var dihedral = findDihedrals(coords).get(0);

		// check a few different angles in [-180,180)
		assertAngle(dihedral, 0.0);
		assertAngle(dihedral, 55.0);
		assertAngle(dihedral, -45.0);
		assertAngle(dihedral, 177.5);
	}

	public static List<DihedralAngle> findDihedrals(AssignedCoords assignedCoords) {
		return assignedCoords.dofs.stream()
			.filter(dof -> dof instanceof DihedralAngle.Dof)
			.map(dof -> ((DihedralAngle.Dof)dof).dihedral)
			.collect(Collectors.toList());
	}

	public static void assertAngle(DihedralAngle dihedral, double angleDegrees) {
		dihedral.setAngle(Math.toRadians(angleDegrees));
		assertThat(dihedral.measureAngleDegrees(), isAbsolutely(angleDegrees, 1e-8));
		assertThat(measureDihedralDegrees(dihedral), isAbsolutely(angleDegrees, 1e-8));
	}

	public static double measureDihedralDegrees(CoordsList coords, int ai, int bi, int ci, int di) {

		var a = new Vector3d();
		var b = new Vector3d();
		var c = new Vector3d();
		var d = new Vector3d();

		coords.get(ai, a);
		coords.get(bi, b);
		coords.get(ci, c);
		coords.get(di, d);

		return Protractor.measureDihedral(
			new double[] { a.x, a.y, a.z }, 0,
			new double[] { b.x, b.y, b.z }, 0,
			new double[] { c.x, c.y, c.z }, 0,
			new double[] { d.x, d.y, d.z }, 0
		);
	}

	public static double measureDihedralDegrees(DihedralAngle dihedral) {
		return measureDihedralDegrees(dihedral.coords.coords, dihedral.ai, dihedral.bi, dihedral.ci, dihedral.di);
	}

	public static double measureDihedralDegrees(CoordsList coords, DihedralAngle dihedral, int di) {
		return measureDihedralDegrees(coords, dihedral.ai, dihedral.bi, dihedral.ci, di);
	}
}
