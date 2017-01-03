package edu.duke.cs.osprey.tools;

import static edu.duke.cs.osprey.TestBase.*;
import static org.junit.Assert.*;

import java.util.Random;

import org.junit.Test;


public class TestDihedrals {
	
	private static final double EpsilonDegrees = 1e-8;
	private static final double EpsilonTrig = 1e-10;
	
	@Test
	public void protractorSquare0() {
		
		double[][] coords = {
			{ 0, 1, 0 },
			{ 0, 0, 0 },
			{ 1, 0, 0 },
			{ 1, 1, 0 }
		};
		
		checkAngle(coords, 0);
	}
	
	@Test
	public void protractorSquare180() {
		
		double[][] coords = {
			{ 0, 1, 0 },
			{ 0, 0, 0 },
			{ 1, 0, 0 },
			{ 1, -1, 0 }
		};
	
		checkAngle(coords, 180);
	}
	
	@Test
	public void protractorSquare90() {
		
		double[][] coords = {
			{ 0, 1, 0 },
			{ 0, 0, 0 },
			{ 1, 0, 0 },
			{ 1, 0, 1 }
		};
		
		checkAngle(coords, 90);
	}
	
	@Test
	public void protractorSquarem90() {
		
		double[][] coords = {
			{ 0, 1, 0 },
			{ 0, 0, 0 },
			{ 1, 0, 0 },
			{ 1, 0, -1 }
		};
		
		checkAngle(coords, -90);
	}
	
	@Test
	public void protractorSquareFineCircle() {
		
		int n = 360*2*1024;
		for (int i=0; i<n; i++) {
			
			// make a square 4-point chain
			double[][] coords = {
				{ 0, 1, 0 },
				{ 0, 0, 0 },
				{ 1, 0, 0 },
				{ 1, 1, 0 }
			};
			checkAngle(coords, 0);
			
			double angleDegrees = (double)i*720/n - 360;
			RotationMatrix rot = new RotationMatrix(1, 0, 0, angleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, angleDegrees);
		}
	}
	
	@Test
	public void protractorFlatlikeFineCircle() {
		
		int n = 360*2*1024;
		for (int i=0; i<n; i++) {
			
			// make a flatlike 4-point chain
			double[][] coords = {
				{ 0, 0.00001, 0 },
				{ 0, 0,     0 },
				{ 1, 0,     0 },
				{ 1, 0.00001, 0 }
			};
			checkAngle(coords, 0);
			
			double angleDegrees = (double)i*720/n - 360;
			RotationMatrix rot = new RotationMatrix(1, 0, 0, angleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, angleDegrees);
		}
	}
	
	@Test
	public void protractorShortAxisFineCircle() {
		
		int n = 360*2*1024;
		for (int i=0; i<n; i++) {
			
			// make a short-axis 4-point chain
			double[][] coords = {
				{ 0,       1, 0 },
				{ 0,       0, 0 },
				{ 0.00001, 0, 0 },
				{ 0.00001, 1, 0 }
			};
			checkAngle(coords, 0);
			
			double angleDegrees = (double)i*720/n - 360;
			RotationMatrix rot = new RotationMatrix(1, 0, 0, angleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, angleDegrees);
		}
	}
	
	@Test
	public void protractorIttyBittyFineCircle() {
		
		int n = 360*2*1024;
		for (int i=0; i<n; i++) {
			
			// make an itty-bitty 4-point chain
			double[][] coords = {
				{ -0.00001,       0.00001, 0 },
				{ 0,        0,       0 },
				{ 0.00001, 0,       0 },
				{ 0.00002, 0.00001, 0 }
			};
			checkAngle(coords, 0);
			
			double angleDegrees = (double)i*720/n - 360;
			RotationMatrix rot = new RotationMatrix(1, 0, 0, angleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, angleDegrees);
		}
	}

	@Test
	public void protractorSkewFineCircle() {
		
		int n = 360*2*1024;
		for (int i=0; i<n; i++) {
			
			double angleDegrees = (double)i*720/n - 360;
			
			// make a skewed 4-point chain
			double[][] coords = {
				{ -1, 1, 0 },
				{ 0, 0, 0 },
				{ 1, 0, 0 },
				{ 2, 1, 0 }
			};
			checkAngle(coords, 0);
			
			RotationMatrix rot = new RotationMatrix(1, 0, 0, angleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, angleDegrees);
		}
	}
	
	@Test
	public void protractorSkewRotFineCircle() {
		
		int n = 360*2*1024;
		for (int i=0; i<n; i++) {
			
			// make a skewed 4-point chain
			double[][] coords = {
				{ -1, 1, 0 },
				{ 0, 0, 0 },
				{ 1, 0, 0 },
				{ 2, 1, 0 }
			};
			checkAngle(coords, 0);
			
			// rotate coords by an arbitrary rotation. the less axis-aligned, the better
			RotationMatrix allRot = new RotationMatrix(3, 8, 1, 36, false);
			for (int j=0; j<4; j++) {
				allRot.applyRotation(coords[j], 0);
			}
			checkAngle(coords, 0);
			
			// then do the dihedral rotation
			double angleDegrees = (double)i*720/n - 360;
			RotationMatrix rot = new RotationMatrix(coords[2][0], coords[2][1], coords[2][2], angleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, angleDegrees);
		}
	}
	
	@Test
	public void protractorSkewRotArbitraryStartFineCircle() {
		
		int n = 360*2*1024;
		for (int i=0; i<n; i++) {
			
			// make a skewed 4-point chain
			double[][] coords = {
				{ -1, 1, 0 },
				{ 0, 0, 0 },
				{ 1, 0, 0 },
				{ 2, 1, 0 }
			};
			checkAngle(coords, 0);
			
			// set the intial dihedral to something arbitrary
			double initialAngleDegrees = -87;
			RotationMatrix rot = new RotationMatrix(coords[2][0], coords[2][1], coords[2][2], initialAngleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, initialAngleDegrees);
			
			// rotate coords by an arbitrary rotation. the less axis-aligned, the better
			RotationMatrix allRot = new RotationMatrix(3, 8, 1, 36, false);
			for (int j=0; j<4; j++) {
				allRot.applyRotation(coords[j], 0);
			}
			
			checkAngle(coords, initialAngleDegrees);
			
			// then do the dihedral rotation
			double nextAngleDegrees = (double)i*720/n - 360;
			double rotAngleDegrees = nextAngleDegrees - Protractor.measureDihedral(coords);
			rot = new RotationMatrix(coords[2][0], coords[2][1], coords[2][2], rotAngleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, nextAngleDegrees);
		}
	}
	
	@Test
	public void protractorSkewRotFineCircleChain() {
		
		// make a skewed 4-point chain
		double[][] coords = {
			{ -1, 1, 0 },
			{ 0, 0, 0 },
			{ 1, 0, 0 },
			{ 2, 1, 0 }
		};
		checkAngle(coords, 0);
		
		// rotate coords by an arbitrary rotation. the less axis-aligned, the better
		RotationMatrix allRot = new RotationMatrix(3, 8, 1, 36, false);
		for (int j=0; j<4; j++) {
			allRot.applyRotation(coords[j], 0);
		}
		checkAngle(coords, 0);
		
		int n = 360*2*1024;
		for (int i=0; i<n; i++) {
			
			// then do the dihedral rotation
			double nextAngleDegrees = (double)i*720/n - 360;
			double rotAngleDegrees = nextAngleDegrees - Protractor.measureDihedral(coords);
			RotationMatrix rot = new RotationMatrix(coords[2][0], coords[2][1], coords[2][2], rotAngleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, nextAngleDegrees);
		}
	}
	
	@Test
	public void protractorSquareRandomRotations() {
		
		// make a square 4-point chain
		double[][] coords = {
			{ 0, 1, 0 },
			{ 0, 0, 0 },
			{ 1, 0, 0 },
			{ 1, 1, 0 }
		};
		checkAngle(coords, 0);
		
		Random rand = new Random(12345);
		
		int n = 1000*1000;
		for (int i=0; i<n; i++) {
			
			// apply a random dihedral rotation
			double nextAngleDegrees = rand.nextDouble()*720 - 360;
			double rotAngleDegrees = nextAngleDegrees - Protractor.measureDihedral(coords);
			RotationMatrix rot = new RotationMatrix(1, 0, 0, rotAngleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, nextAngleDegrees);
		}
	}
	
	@Test
	public void reallyPileOnTheRoundoffError() {
		
		// make an itty-bitty 4-point chain
		double[][] coords = {
			{ -0.00001, 0.00001, 0 },
			{  0,       0,       0 },
			{  0.00001, 0,       0 },
			{  0.00002, 0.00001, 0 }
		};
		checkAngle(coords, 0);
		
		// rotate coords by an arbitrary rotation. the less axis-aligned, the better
		RotationMatrix allRot = new RotationMatrix(3, 8, 1, -83, false);
		for (int j=0; j<4; j++) {
			allRot.applyRotation(coords[j], 0);
		}
		checkAngle(coords, 0);
		
		Random rand = new Random(12345);
		
		int n = 360*2*1024;
		for (int i=0; i<n; i++) {
			
			// apply a random dihedral rotation
			RotationMatrix rot = new RotationMatrix(coords[2][0], coords[2][1], coords[2][2], rand.nextDouble()*360, false);
			rot.applyRotation(coords[3], 0);
			
			// then do the dihedral rotation
			double nextAngleDegrees = (double)i*720/n - 360;
			double rotAngleDegrees = nextAngleDegrees - Protractor.measureDihedral(coords);
			rot = new RotationMatrix(coords[2][0], coords[2][1], coords[2][2], rotAngleDegrees, false);
			rot.applyRotation(coords[3], 0);
			
			checkAngle(coords, nextAngleDegrees);
		}
	}
	
	private void checkAngle(double[][] coords, double angleDegrees) {
		checkAngle(coords, angleDegrees, EpsilonDegrees, EpsilonTrig);
	}
	
	private void checkAngle(double[][] coords, double angleDegrees, double epsilonDegrees, double epsilonTrig) {
		angleDegrees = Protractor.normalizeDegrees(angleDegrees);
		assertThat(Protractor.measureDihedral(coords), isDegrees(angleDegrees, epsilonDegrees));
		double angleRadians = Math.toRadians(angleDegrees);
		assertThat(Protractor.measureDihedralSinCos(coords), isAbsolutely(new double[] { Math.sin(angleRadians), Math.cos(angleRadians) }, epsilonTrig));
	}
	
	//@Test
	// NOTE: this test always fails because half angle computations are not numerically stable!
	// don't use them in practice!
	public void halfAngles() {
		
		int n = 360*2;
		for (int i=0; i<n; i++) {
			
			double angleDegrees = (double)i*720/n - 360;
			double s = Math.sin(Math.toRadians(angleDegrees));
			double c = Math.cos(Math.toRadians(angleDegrees));
			
			double halfAngleDegrees = Protractor.normalizeDegrees(angleDegrees)/2;
			double hs = Math.sin(Math.toRadians(halfAngleDegrees));
			double hc = Math.cos(Math.toRadians(halfAngleDegrees));
			
			@SuppressWarnings("deprecation")
			double[] halfSinCos = Protractor.getHalfAngleSinCos(s, c);
			
			/* DEBUG
			System.out.println(String.format("%7.2f %7.2f %16.12f %2d %16.12f %2d   /2   %7.2f %16.12f %16.12f   ->   %16.12f   %16.12f",
				angleDegrees, Protractor.normalizeDegrees(angleDegrees), c, (int)Math.signum(c), s, (int)Math.signum(s),
				halfAngleDegrees, hc, hs,
				halfSinCos[1], halfSinCos[0]
			));
			*/
			
			assertThat(halfSinCos[0], isAbsolutely(hs, EpsilonDegrees));
			assertThat(halfSinCos[1], isAbsolutely(hc, EpsilonDegrees));
		}
	}
}
