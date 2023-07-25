package edu.duke.cs.osprey.kstar.compiled;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.MatcherAssert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.ematrix.compiled.ErefCalculator;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

import java.util.HashMap;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Most just regression tests to make sure a bug don't cause lower energy bounds on conformations to become looser.
 * Those errors are harder to catch because they don't cause K* calculations to fail, just take MUCH MUCH longer.
 * So these tests turn bounds changes into explicit errors, but they really only work as regression tests.
 */
public class TestBounds {

	private boolean generate = false;

	@Test
	public void test2RL0Complex() {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.complex.ccsx"));
		var posInterDist = PosInterDist.DesmetEtAl1992;

		var bounds = new HashMap<String,Double>();
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.796);
		bounds.put("156 PHE=ALA 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.802);
		bounds.put("156 PHE=ILE 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.814);
		bounds.put("156 PHE=LEU 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.823);
		bounds.put("156 PHE=TYR 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.804);
		bounds.put("156 PHE=VAL 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.810);
		bounds.put("156 PHE=phe 172 LYS=ASP 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.096);
		bounds.put("156 PHE=phe 172 LYS=GLU 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.621);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ALA 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.132);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=LEU 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",    21.399);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=PHE 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",    94.842);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=TYR 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",   104.555);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=VAL 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.145);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=ASN 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.422);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=SER 649 PHE=phe 650 ASP=asp 651 GLU=glu",     2.238);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=ALA 650 ASP=asp 651 GLU=glu",     2.138);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=ILE 650 ASP=asp 651 GLU=glu",     2.506);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=LEU 650 ASP=asp 651 GLU=glu",     2.750);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=TYR 650 ASP=asp 651 GLU=glu",     3.395);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=VAL 650 ASP=asp 651 GLU=glu",     2.250);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=GLU 651 GLU=glu",     1.828);
		bounds.put("156 PHE=phe 172 LYS=lys 192 ILE=ile 193 THR=thr 649 PHE=phe 650 ASP=asp 651 GLU=ASP",     2.186);

		checkBounds(complex, posInterDist, false, bounds);
	}

	private void checkBounds(ConfSpace confSpace, PosInterDist posInterDist, boolean includeEref, HashMap<String,Double> bounds) {

		final boolean includeStaticStatic = true;

		var parallelism = Parallelism.makeCpu(4);
		try (var tasks = parallelism.makeTaskExecutor()) {

			//ConfEnergyCalculator ecalc = ConfEnergyCalculator.makeBest(confSpace, parallelism);
			ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(confSpace);

			// calc the reference energies if needed
			SimpleReferenceEnergies eref = null;
			if (includeEref) {
				eref = new ErefCalculator.Builder(ecalc)
					.build()
					.calc(tasks);
			}

			// calc the emat
			EnergyMatrix emat = new EmatCalculator.Builder(ecalc)
				.setPosInterDist(posInterDist)
				.setReferenceEnergies(eref)
				.setIncludeStaticStatic(includeStaticStatic)
				.build()
				.calc(tasks);

			for (var seq : confSpace.seqSpace.getSequences(1)) {

				// use A* to get the lowest-scoring sequence
				var rcs = seq.makeRCs(confSpace);
				var astar = new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build();
				var conf = astar.nextConf();

				// calc the energy
				var posInters = PosInterDist.all(confSpace, eref, conf.getAssignments());
				double energy = ecalc.minimizeEnergy(conf.getAssignments(), posInters);

				double gap = energy - conf.getScore();

				if (generate) {
					log("bounds.put(\"%s\", %9.3f);", seq.toString(), gap);
				} else {
					double gapExp = bounds.get(seq.toString());
					assertThat(gap, isAbsolutely(gapExp, 1e-1));
				}
			}
		}
	}

	public static void main(String[] args) {
		var tester = new TestBounds();
		tester.generate = true;
		tester.test2RL0Complex();
	}
}
