package edu.duke.cs.osprey.kstar.pfunc;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Samples the 2D decision space of estimating a partition function
 * (ie, refine the upper bound, or refine the lower bound)
 * and reports the delta value of the pfunc estimation as a VTK file.
 *
 * Lawrence Livermore National Lab's VisIt software is a good tool for
 * visualizing this data:
 * https://wci.llnl.gov/simulation/computer-codes/visit/
 *
 * This class also tracks "traces" for visualization via VTK files.
 * ie, paths over this surface taken by a partition function calculator.
 */
public class PfuncSurface {

	public final int scoreBatch;
	public final int numScoreBatches;
	public final int numEnergies;
	public final int numScores;

	private double[][] delta = null;
	private List<Trace> traces = new ArrayList<>();

	public PfuncSurface(int scoreBatch, int numScoreBatches, int numEnergies) {
		this.scoreBatch = scoreBatch;
		this.numScoreBatches = numScoreBatches;
		this.numEnergies = numEnergies;
		this.numScores = scoreBatch*numScoreBatches;
	}

	public void sample(ConfEnergyCalculator confEcalc, EnergyMatrix emat) {
		sample(confEcalc, emat, new RCs(confEcalc.confSpace));
	}

	public void sample(ConfEnergyCalculator confEcalc, EnergyMatrix emat, RCs rcs) {

		// collect all the scores in a list
		log("collecting scores...");
		List<ConfSearch.ScoredConf> scoredConfs = new ArrayList<>(numScores);
		{
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build();
			for (int i=0; i<numScoreBatches*scoreBatch; i++) {
				scoredConfs.add(astar.nextConf());
			}
		}

		BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

		// collect all the energies in a list
		log("computing energies...");
		List<ConfSearch.EnergiedConf> energiedConfs = confEcalc.calcAllEnergies(scoredConfs.subList(0, numEnergies));

		// and their Boltzmann-weighted energies
		log("computing Boltzmann weights...");
		List<BigDecimal> scoreWeights = new ArrayList<>(numEnergies);
		for (ConfSearch.ScoredConf conf : scoredConfs) {
			scoreWeights.add(bcalc.calc(conf.getScore()));
		}
		List<BigDecimal> energyWeights = new ArrayList<>(energiedConfs.size());
		for (ConfSearch.EnergiedConf econf : energiedConfs) {
			energyWeights.add(bcalc.calc(econf.getEnergy()));
		}

		class State {

			BigDecimal numConfs = new BigDecimal(rcs.getNumConformations());
			int numScoredConfs = 0;
			BigDecimal upperScoreWeightSum = BigDecimal.ZERO;
			BigDecimal lastSWeight = BigDecimal.ZERO;

			BigDecimal lowerScoreWeightSum = BigDecimal.ZERO;
			BigDecimal energyWeightSum = BigDecimal.ZERO;

			double calcDelta() {
				BigDecimal unscoredBound = numConfs
					.subtract(BigDecimal.valueOf(numScoredConfs))
					.multiply(lastSWeight);
				BigDecimal adjustedUpperBound = upperScoreWeightSum
					.subtract(lowerScoreWeightSum)
					.add(energyWeightSum)
					.add(unscoredBound);
				return MathTools.bigDivide(
					adjustedUpperBound.subtract(energyWeightSum),
					adjustedUpperBound,
					PartitionFunction.decimalPrecision
				).doubleValue();
			}
		}
		State state = new State();

		// initialize delta
		delta = new double[numScoreBatches + 1][numEnergies + 1];
		for (double[] d : delta) {
			Arrays.fill(d, 1.0);
		}

		// sample the pfunc surface
		log("sampling...");
		for (int s=1; s<=numScoreBatches; s++) {

			// go through a score batch and record the delta
			for (int i=0; i<scoreBatch; i++) {
				BigDecimal scoreWeight = scoreWeights.get((s - 1)*scoreBatch + i);
				state.upperScoreWeightSum = state.upperScoreWeightSum.add(scoreWeight);
				state.lastSWeight = scoreWeight;
			}
			state.numScoredConfs += scoreBatch;

			state.lowerScoreWeightSum = BigDecimal.ZERO;
			state.energyWeightSum = BigDecimal.ZERO;

			delta[s][0] = state.calcDelta();

			// go through all the energies and record the deltas
			for (int e=1; e<=Math.min(s*scoreBatch, numEnergies); e++) {

				BigDecimal scoreWeight = scoreWeights.get(e - 1);
				BigDecimal energyWeight = energyWeights.get(e - 1);

				state.lowerScoreWeightSum = state.lowerScoreWeightSum.add(scoreWeight);
				state.energyWeightSum = state.energyWeightSum.add(energyWeight);

				delta[s][e] = state.calcDelta();
			}
		}

		log("surface sampling complete!");
	}

	public void write(File file)
	throws Exception {

		/* write out the vtk file, e.g.:
			# vtk DataFile Version 3.0
			beer is super awesome!
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
		try (Writer out = new FileWriter(file)) {

			// write headers
			out.write("# vtk DataFile Version 3.0\n");
			out.write("whatever\n");
			out.write("ASCII\n");
			out.write("DATASET RECTILINEAR_GRID\n");
			out.write(String.format("DIMENSIONS %d %d %d\n", numScoreBatches + 1, numEnergies + 1, 1));

			out.write(String.format("X_COORDINATES %d float\n", numScoreBatches + 1));
			for (int i=0; i<=numScoreBatches; i++) {
				if (i % 10 > 0) {
					out.write(" ");
				} else if (i > 0) {
					out.write("\n");
				}
				// use confs instead of time
				//out.write(String.format("%.2f", scoreMs[i]));
				out.write(String.format("%d", i));
			}
			out.write("\n");

			out.write(String.format("Y_COORDINATES %d float\n", numEnergies + 1));
			for (int i=0; i<=numEnergies; i++) {
				if (i % 10 > 0) {
					out.write(" ");
				} else if (i > 0) {
					out.write("\n");
				}
				// use confs instead of time
				//out.write(String.format("%.2f", energyMs[i]));
				out.write(String.format("%d", i));
			}
			out.write("\n");

			out.write("Z_COORDINATES 1 float\n");
			out.write("0\n");

			out.write(String.format("POINT_DATA %d\n", (numScoreBatches + 1)*(numEnergies + 1)));
			out.write("FIELD fieldDelta 1\n");
			out.write(String.format("delta 1 %d float\n", (numScoreBatches + 1)*(numEnergies + 1)));
			for (int e=0; e<=numEnergies; e++) {
				for (int s=0; s<=numScoreBatches; s++) {
					if (s % 10 > 0) {
						out.write(" ");
					} else if (s > 0) {
						out.write("\n");
					}
					out.write(String.format("%.4f", delta[s][e]));
				}
				out.write("\n");
			}
		}
		log("wrote to file %s", file.getAbsolutePath());
	}

	public class Trace {

		class Point {
			int scores;
			int energies;
			double delta;
		}

		public final List<Point> points = new ArrayList<>();

		public Trace() {
			traces.add(this);
		}

		public void step(int scores, int energies, double delta) {
			Point p = new Point();
			p.scores = scores;
			p.energies = energies;
			p.delta = delta;
			points.add(p);
		}
	}

	public void writeTraces(File file) {

		/* write the trace as a VTK file, e.g.:
			# vtk DataFile Version 3.0
			beer is super awesome!
			ASCII
			DATASET POLYDATA
			POINTS 3 float
			0 0 0
			1 1 1
			2 2 2
			LINES 2 6
			2 0 1
			2 1 2
		*/
		try (Writer out = new FileWriter(file)) {

			// write headers
			out.write("# vtk DataFile Version 3.0\n");
			out.write("whatever\n");
			out.write("ASCII\n");
			out.write("DATASET POLYDATA\n");

			int numPoints = traces.stream().mapToInt((trace) -> trace.points.size()).sum();
			out.write(String.format("POINTS %d float\n", numPoints));
			for (Trace trace : traces) {
				for (Trace.Point p : trace.points) {
					out.write(String.format("%d %d %d\n", p.scores/scoreBatch, p.energies, 0));
				}
			}

			int numLines = numPoints - traces.size();
			out.write(String.format("LINES %d %d\n", numLines, numLines*3));
			int offset = 0;
			for (Trace trace : traces) {
				for (int i=0; i<trace.points.size()-1; i++) {
					out.write(String.format("2 %d %d\n", offset + i, offset + i + 1));
				}
				offset += trace.points.size();
			}

			out.write(String.format("POINT_DATA %d\n", numPoints));
			out.write("FIELD fieldDelta 1\n");
			out.write(String.format("delta 1 %d float\n", numPoints));
			for (Trace trace : traces) {
				for (int i=0; i<trace.points.size(); i++) {
					if (i % 10 > 0) {
						out.write(" ");
					} else if (i > 0) {
						out.write("\n");
					}
					out.write(String.format("%.4f", trace.points.get(i).delta));
				}
			}
			out.write("\n");

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}
}
