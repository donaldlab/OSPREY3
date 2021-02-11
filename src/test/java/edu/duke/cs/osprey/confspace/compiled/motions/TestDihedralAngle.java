package edu.duke.cs.osprey.confspace.compiled.motions;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.isAbsolutely;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Protractor;
import org.joml.Vector3d;
import org.junit.Test;


public class TestDihedralAngle {

	@Test
	public void test() {

		// load a conf space
		var cs = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/dipeptide.5hydrophobic.ccsx"));

		// make a VAL-GLY conformation
		var conf = new int[] { 14, 1 };
		assertThat(cs.confType(0, conf[0]), is("VAL"));
		assertThat(cs.confType(1, conf[1]), is("GLY"));
		var coords = cs.makeCoords(conf);

		// VAL should have just one dihedral angle
		assertThat(coords.dofs.size(), is(1));
		var dof = (DihedralAngle.Dof)coords.dofs.get(0);
		var dihedral = dof.dihedral;

		// check a few different angles in [-180,180)
		assertAngle(dihedral, 0.0);
		assertAngle(dihedral, 55.0);
		assertAngle(dihedral, -45.0);
		assertAngle(dihedral, 177.5);
	}

	private static void assertAngle(DihedralAngle dihedral, double angleDegrees) {
		dihedral.setAngle(Math.toRadians(angleDegrees));
		assertThat(dihedral.measureAngleDegrees(), isAbsolutely(angleDegrees, 1e-8));
		assertThat(measureDihedralDegrees(dihedral), isAbsolutely(angleDegrees, 1e-8));
	}

	private static double measureDihedralDegrees(DihedralAngle dihedral) {

		var a = new Vector3d();
		var b = new Vector3d();
		var c = new Vector3d();
		var d = new Vector3d();

		dihedral.coords.coords.get(dihedral.ai, a);
		dihedral.coords.coords.get(dihedral.bi, b);
		dihedral.coords.coords.get(dihedral.ci, c);
		dihedral.coords.coords.get(dihedral.di, d);

		return Protractor.measureDihedral(
			new double[] { a.x, a.y, a.z }, 0,
			new double[] { b.x, b.y, b.z }, 0,
			new double[] { c.x, c.y, c.z }, 0,
			new double[] { d.x, d.y, d.z }, 0
		);
	}
}
