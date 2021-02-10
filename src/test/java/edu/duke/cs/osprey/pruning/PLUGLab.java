/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.pruning;


import edu.duke.cs.osprey.Benchmark;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.*;
import edu.duke.cs.osprey.tools.Profiler;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


public class PLUGLab {

	public static void main(String[] args) {

		Parallelism parallelism = Parallelism.makeCpu(4);

		// load a protein
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A23").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // asn
		strand.flexibility.get("A27").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // thr
		strand.flexibility.get("A36").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // ile
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// benchmark PLUG
		RCTuple tuple = new RCTuple(1, 0);
		PLUG plug = new PLUG(confSpace);
		Benchmark b = new Benchmark(10, 200, () -> {
			plug.shouldPruneTuple(tuple, 0.4);
		});
		log("tuple %s: %s", tuple, b);

		// TEMP
		if (true) return;

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(parallelism)
			.build()) {

			EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(ecalc)
				.setIsMinimizing(false)
				.build();

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				//.setEnergyPartition(EnergyPartition.AllOnPairs)
				.build();

			// TEMP: look at conf [7, 12, 7]
			// looks like minimization doesn't avoid the clash, but still gets a relatively good energy
			// clashes:
			//    A27:HG23 <-> A30:OE1               overlap=   0.408  >0.000 BadClash
			//    A27:C    <-> A27:HG22              overlap=   0.551  >0.000 BadClash
			int[] assignments = new int[] { 7, 12, 7 };
			RCTuple t = new RCTuple(assignments);
			ResidueInteractions inters = ResInterGen.of(confSpace)
				.addIntras(t)
				.addInters(t)
				.addShell(t)
				.make();

			EnergyCalculator.EnergiedParametricMolecule rigidEpmol = rigidEcalc.calcEnergy(confSpace.makeMolecule(t), inters);
			log("before: %.4f", rigidEpmol.energy);
			EnergyCalculator.EnergiedParametricMolecule minEpmol = ecalc.calcEnergy(confSpace.makeMolecule(t), inters);
			log("after:  %.4f", minEpmol.energy);

			Probe probe = new Probe();
			probe.matchTemplates(confSpace);
			probe.getInteractions(minEpmol.pmol.mol.residues, inters, ecalc.resPairCache.connectivity).stream()
				.filter(interaction -> interaction.contact.isContact)
				.sorted(Comparator.comparing(interaction -> interaction.getViolation(0.0)))
				.forEach(interaction -> log("%s   %s", interaction.atomPair, interaction));

			log("%s", new ResidueForcefieldBreakdown.ByPosition(confEcalc, assignments, minEpmol).breakdownForcefield(ResidueForcefieldBreakdown.Type.VanDerWaals));
			log("%d -> %s", 25, minEpmol.pmol.mol.residues.get(25).getPDBResNumber());
			log("%s", new ResidueForcefieldBreakdown.ByResidue(confEcalc, minEpmol).breakdownForcefield(ResidueForcefieldBreakdown.Type.VanDerWaals));

			PDBIO.writeFile(minEpmol, new File("clash.pdb"));

			if (true) return;

			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build()
				.calcEnergyMatrix();

			// get the best few confs
			List<ConfSearch.ScoredConf> confs = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.build()
				.nextConfs(100);
			List<ConfSearch.EnergiedConf> econfs = confEcalc.calcAllEnergies(confs);
			econfs.sort(Comparator.comparing((econf) -> econf.getEnergy()));
			log("before pruning:");
			for (int i=0; i<econfs.size(); i++) {
				log("%s %.4f", Conf.toString(econfs.get(i).getAssignments()), econfs.get(i).getEnergy());
			}

			// run PLUG
			Stopwatch pmatStopwatch = new Stopwatch().start();
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setParallelism(parallelism)
				.setThreshold(null) // no steric pruning
				.setSinglesPlugThreshold(0.4)
				.setPairsPlugThreshold(0.4)
				.setShowProgress(true)
				.run(confSpace, null);
			log("%s", pmatStopwatch.stop().getTime(2));

			// get the best few confs again
			confs = new ConfAStarTree.Builder(emat, pmat)
				.setTraditional()
				.build()
				.nextConfs(100);
			econfs = confEcalc.calcAllEnergies(confs);
			econfs.sort(Comparator.comparing((econf) -> econf.getEnergy()));
			log("after pruning:");
			for (int i=0; i<econfs.size(); i++) {
				log("%s %.4f", Conf.toString(econfs.get(i).getAssignments()), econfs.get(i).getEnergy());
			}
		}
	}

	/** try to find a boundary by brute force sampling */
	private static boolean containsBoundary(ParametricMolecule pmol, Probe.AtomPair pair, double tolerance, int numSamplesPerDof) {

		// make voxel properties easy to access
		int n = pmol.dofs.size();
		double[] vmin = new double[n];
		double[] vmax = new double[n];
		for (int d=0; d<n; d++) {
			vmin[d] = pmol.dofBounds.getMin(d);
			vmax[d] = pmol.dofBounds.getMax(d);
		}

		boolean foundNeg = false;
		boolean foundPos = false;

		double[] i = new double[n];
		Arrays.fill(i, 0);

		while (true) {

			for (int d=0; d<n; d++) {
				double xd = vmin[d] + (vmax[d] - vmin[d])*i[d]/(numSamplesPerDof - 1);
				pmol.dofs.get(d).apply(xd);
			}

			double violation = pair.getViolation(tolerance);
			if (violation < 0) {
				foundNeg = true;
			} else if (violation > 0) {
				foundPos = true;
			} else { // violation == 0
				return true;
			}

			// advance to next indices
			i[0]++;
			for (int d=0; d<n; d++) {
				if (i[d] >= numSamplesPerDof) {

					if (d == n - 1) {
						// we're done
						return foundNeg && foundPos;
					}

					i[d + 1]++;

					for (int d2=0; d2<=d; d2++) {
						i[d] = 0;
					}
				}
			}
		}
	}

	/**
	 * sample two degrees of freedom for a parametric molecule and write the violations for an atom pair to a file
	 *
	 * Lawrence Livermore National Lab's VisIt software is a good tool for
	 * visualizing this data:
	 * https://wci.llnl.gov/simulation/computer-codes/visit/
	 *
	 * (use a pseudocolor plot with a contour at 0)
	 */
	public static void writeViolations2D(ParametricMolecule pmol, int d1, int d2, Probe.AtomPair pair, double tolerance, File file) {

		// figure out the dofs
		DegreeOfFreedom dof1 = pmol.dofs.get(d1);
		double min1 = pmol.dofBounds.getMin(d1);
		double max1 = pmol.dofBounds.getMax(d1);
		int samples1 = 40;
		Function<Integer,Double> sample1 = (i1) -> min1 + (max1 - min1)*i1/(samples1 - 1);

		DegreeOfFreedom dof2 = pmol.dofs.get(d2);
		double min2 = pmol.dofBounds.getMin(d2);
		double max2 = pmol.dofBounds.getMax(d2);
		int samples2 = 40;
		Function<Integer,Double> sample2 = (i2) -> min2 + (max2 - min2)*i2/(samples2 - 1);

		try (Writer out = new FileWriter(file)) {

			/* write out the vtk file, e.g.:
				# vtk DataFile Version 3.0
				any text here
				ASCII
				DATASET RECTILINEAR_GRID
				DIMENSIONS 3 4 1
				X_COORDINATES 3 float
				0 2 4
				Y_COORDINATES 4 float
				1 2 3 4
				Z_COORDINATES 1 float
				0
				POINT_DATA 12
				FIELD FieldData 1
				cellscalar 1 12 float
				0 0 0
				1 1 1
				2 2 2
				3 3 3

				see: http://www.visitusers.org/index.php?title=ASCII_VTK_Files
				and: https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
			*/

			// write headers
			out.write("# vtk DataFile Version 3.0\n");
			out.write("whatever\n");
			out.write("ASCII\n");
			out.write("DATASET RECTILINEAR_GRID\n");
			out.write(String.format("DIMENSIONS %d %d %d\n", samples1, samples2, 1));

			// write x axis (x1)
			out.write(String.format("X_COORDINATES %d float\n", samples1));
			for (int i1=0; i1<samples1; i1++) {
				if (i1 % 10 > 0) {
					out.write(" ");
				} else if (i1 > 0) {
					out.write("\n");
				}
				out.write(String.format("%.2f", sample1.apply(i1)));
			}
			out.write("\n");

			// write y axis (x2)
			out.write(String.format("Y_COORDINATES %d float\n", samples2));
			for (int i2=0; i2<samples2; i2++) {
				if (i2 % 10 > 0) {
					out.write(" ");
				} else if (i2 > 0) {
					out.write("\n");
				}
				out.write(String.format("%.2f", sample2.apply(i2)));
			}
			out.write("\n");

			out.write("Z_COORDINATES 1 float\n");
			out.write("0\n");

			out.write(String.format("POINT_DATA %d\n", samples1*samples2));
			out.write("FIELD fieldDelta 1\n");
			out.write(String.format("delta 1 %d float\n", samples1*samples2));

			// sample DoFs (in reverse axis order, for the vtk file)
			int p = 0;
			for (int i2=0; i2<samples2; i2++) {
				double x2 = sample2.apply(i2);
				dof2.apply(x2);

					for (int i1=0; i1<samples1; i1++) {
					double x1 = sample1.apply(i1);
					dof1.apply(x1);

					// calculate the violation
					double violation = pair.getViolation(tolerance);

					// write the point
					if (p % 10 > 0) {
						out.write(" ");
					} else if (p > 0) {
						out.write("\n");
					}
					p++;
					out.write(String.format("%.4f", violation));
				}
				out.write("\n");
			}

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}
}
