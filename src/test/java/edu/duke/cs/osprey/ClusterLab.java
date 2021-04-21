package edu.duke.cs.osprey;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.ForkCluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


public class ClusterLab {

	public static void main(String[] args) {

		Cluster.fixHazelcastLogging();

		//forkCluster();
		//multiProcessCluster(args);
		slurmCluster();
	}

	private static void forkCluster() {
		ForkCluster.run(2, true, ClusterLab.class, ClusterLab::run);
	}

	private static void multiProcessCluster(String[] args) {

		// parse the args to get nodeId and cluster size
		String jobId = args[0];
		int nodeId = Integer.parseInt(args[1]);
		int numNodes = Integer.parseInt(args[2]);
		boolean clientIsMember = false;

		run(new Cluster("MultiProcessCluster", jobId, nodeId, numNodes, clientIsMember));
	}

	private static void slurmCluster() {
		run(Cluster.fromSLURM(true));
	}

	private static void run(Cluster cluster) {

		Stopwatch stopwatch = new Stopwatch().start();

		Parallelism parallelism = Parallelism.makeCpu(4);

		// set up a toy design
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.asList(), confSpaces.ffparams)
			.setCluster(cluster)
			.setParallelism(parallelism)
			.build()) {

			// run K*
			KStarScoreWriter.Formatter testFormatter = (KStarScoreWriter.ScoreInfo info) -> {

				Function<PartitionFunction.Result,String> formatPfunc = (pfuncResult) -> {
					if (pfuncResult.status == PartitionFunction.Status.Estimated) {
						return String.format("%12e", pfuncResult.values.qstar.doubleValue());
					}
					return "null";
				};

				return String.format("assertSequence(result, %3d, \"%s\", %-12s, %-12s, %-12s, epsilon); // protein %s ligand %s complex %s K* = %s",
					info.sequenceNumber,
					info.sequence.toString(Sequence.Renderer.ResType),
					formatPfunc.apply(info.kstarScore.protein),
					formatPfunc.apply(info.kstarScore.ligand),
					formatPfunc.apply(info.kstarScore.complex),
					info.kstarScore.protein.toString(),
					info.kstarScore.ligand.toString(),
					info.kstarScore.complex.toString(),
					info.kstarScore.toString()
				);
			};

			// configure K*
			KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(0.95)
				.setStabilityThreshold(null)
				.addScoreConsoleWriter(testFormatter)
				.setMaxSimultaneousMutations(1)
				.build();
			KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				// turn off the default confdb for tests
				info.confDBFile = null;

				SimpleConfSpace confSpace = (SimpleConfSpace)info.confSpace;

				// compute reference energies locally on each node
				SimpleReferenceEnergies eref;
				try (EnergyCalculator localEcalc = ecalc.local()) {
					eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, localEcalc)
						.build()
						.calcReferenceEnergies();
				}

				// how should we define energies of conformations?
				info.confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.setReferenceEnergies(eref)
					.build();

				// calc the energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.build()
					.calcEnergyMatrix();

				info.pfuncFactory = rcs -> new GradientDescentPfunc(
					info.confEcalc,
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					rcs.getNumConformations()
				);
			}

			// run K*
			kstar.run(ecalc.tasks);
		}

		if (cluster.nodeId > 0) {
			log("MEMBER %d: finished in %s", cluster.nodeId, stopwatch.stop().getTime(2));
		} else {
			log("CLIENT: finished in %s", stopwatch.stop().getTime(2));
		}
	}
}
