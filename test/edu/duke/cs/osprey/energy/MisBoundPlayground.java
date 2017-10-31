package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Progress;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class MisBoundPlayground {

	public static void main(String[] args)
	throws Exception {

		ForcefieldParams ffparams = new ForcefieldParams();

		// load non-natural AA stuff
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.addTemplates(FileTools.readFile("all_nuc94_and_gr.in"))
			.addTemplateCoords(FileTools.readFile("all_amino_coords.in"))
			.addRotamers(FileTools.readFile("GenericRotamers.dat"))
			.build();

		Strand strand = new Strand.Builder(PDBIO.readFile("model0_antibody.pdb"))
			.setTemplateLibrary(templateLib)
			.build();

		List<String> mutations = Arrays.asList(
			Strand.WildType,
			"ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "CYX", "MET", "SER",
			"THR", "LYS", "ARG", "HIP", "HIE", "ASP", "GLU", "ASN", "GLN", "GLY"
		);

		for (String resNum : Arrays.asList("H1901", "H1904", "H1905")) {
			strand.flexibility.get(resNum).setLibraryRotamers(mutations).addWildTypeRotamers().setContinuous();
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.setShellDistance(6)
			.build();

		// which bound to use?
		//EnergyPartition epart = EnergyPartition.Traditional;
		EnergyPartition epart = EnergyPartition.AllOnPairs;

		new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.make(4, 1, 8))
			//.setType(EnergyCalculator.Type.Cpu)
			.setType(EnergyCalculator.Type.ResidueCudaCCD)
			.use((ecalc) -> {

				SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcReferenceEnergies();

				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.setReferenceEnergies(eref)
					.setEnergyPartition(epart)
					.build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
					.setCacheFile(new File("emat.misbound.shell." + epart.name().toLowerCase() + ".dat"))
					.build()
					.calcEnergyMatrix();

				ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
					.setTraditional()
					.setShowProgress(true)
					.build();

				// ok, all the pre-reqs are setup, what to do now?
				//doSearch(astar, confEcalc);
				//analyzeMisbounds(confSpace, confEcalc, emat);
				findMisbounds(confSpace, astar, confEcalc, emat);
			});
	}

	public static void doSearch(ConfAStarTree astar, ConfEnergyCalculator confEcalc) {
		new SimpleGMECFinder.Builder(astar, confEcalc)
			.build()
			.find(5.0);
	}

	public static void analyzeMisbounds(SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc, EnergyMatrix emat) {

		List<int[]> misboundedConfsTraditional = Arrays.asList(
			new int[] { 22, 131, 159 }, // Conformation score (-14.511145) is not a lower bound on the energy (-14.869177).
			new int[] { 27, 131, 159 } // Conformation score (-13.803699) is not a lower bound on the energy (-14.221509).
		);

		List<int[]> misboundedConfsAllOnPairs = Arrays.asList(
			new int[] { 22, 140, 149 }, // Conformation score (-14.033345) is not a lower bound on the energy (-14.211877).
			new int[] { 22, 68, 149 }, // Conformation score (-13.812986) is not a lower bound on the energy (-13.971702).
			new int[] { 22, 129, 149 }, // Conformation score (-13.759669) is not a lower bound on the energy (-13.866282).
			new int[] { 22, 70, 149 }, // Conformation score (-13.691940) is not a lower bound on the energy (-13.805368).
			new int[] { 22, 131, 177 }, // Conformation score (-13.665316) is not a lower bound on the energy (-13.847284).
			new int[] { 22, 74, 149 } // Conformation score (-13.643942) is not a lower bound on the energy (-13.784523).
		);

		// check this conf
		int[] conf = { 22, 131, 159 };

		// break down the forcefield and bound
		ResidueForcefieldBreakdown.ByPosition breakdown = new ResidueForcefieldBreakdown.ByPosition(confEcalc, conf);
		EnergyMatrix forcefieldBreakdown = breakdown.breakdownForcefield(ResidueForcefieldBreakdown.Type.All);
		EnergyMatrix scoreBreakdown = breakdown.breakdownScore(emat);

		System.out.println(String.format("Mis-bounded conf:"
				+ "\n\tRCs:       %s"
				+ "\n\tSequence:  %s"
				+ "\n\tRotamers:  %s"
				+ "\n\tScore:     %12.6f"
				+ "\n\tEnergy:    %12.6f",
			confSpace.formatConfRCs(conf),
			confSpace.formatConfSequence(conf),
			confSpace.formatConfRotamers(conf),
			scoreBreakdown.sum(),
			forcefieldBreakdown.sum()
		));
		System.out.println("forcefield: " + forcefieldBreakdown);
		System.out.println("score: " + scoreBreakdown);
	}

	public static void findMisbounds(SimpleConfSpace confSpace, ConfAStarTree astar, ConfEnergyCalculator confEcalc, EnergyMatrix emat)
	throws Exception {

		// get a bunch of confs to check
		List<ConfSearch.ScoredConf> confs = new ArrayList<>();
		//int numConfs = 100; // for testing
		int numConfs = 18000; // for Traditional
		//int numConfs = 12800; // for AllOnPairs
		for (int i=0; i<numConfs; i++) {
			ConfSearch.ScoredConf conf = astar.nextConf();
			if (conf == null) {
				break;
			}
			confs.add(conf);
		}

		/*
		List<ConfSearch.ScoredConf> confs = Arrays.asList(
			new ConfSearch.ScoredConf(new int[] { 22, 131, 159 }, 0.0)
		);
		*/

		double errorThreshold = 1e-6;
		File outFile = new File("misbounds." + confEcalc.epart.name().toLowerCase() + ".txt");

		// clear the out file
		outFile.delete();

		// check the confs for mis-bounds
		Progress progress = new Progress(numConfs);
		for (ConfSearch.ScoredConf conf : confs) {

			confEcalc.tasks.submit(() -> {

				// break down the forcefield and bound
				ResidueForcefieldBreakdown.ByPosition breakdown = new ResidueForcefieldBreakdown.ByPosition(confEcalc, conf.getAssignments());
				EnergyMatrix forcefieldBreakdown = breakdown.breakdownForcefield(ResidueForcefieldBreakdown.Type.All);
				EnergyMatrix scoreBreakdown = breakdown.breakdownScore(emat);

				// are any entries mis-bounded?
				List<String> reports = new ArrayList<>();
				for (SimpleConfSpace.Position pos1 : confSpace.positions) {

					{
						double boundEnergy = scoreBreakdown.getOneBody(pos1.index, 0);
						double forcefieldEnergy = forcefieldBreakdown.getOneBody(pos1.index, 0);
						if (boundEnergy > forcefieldEnergy + errorThreshold) {
							reports.add(String.format("%-30s    cell-score %12.6f    Sequence: %s    cell-energy %12.6f    gap: %12.6f    conf-score: %12.6f    conf-energy: %12.6f    gap: %12.6f",
								"Single: " + pos1.formatConfPos(conf),
								boundEnergy,
								confSpace.formatConf(conf),
								forcefieldEnergy,
								boundEnergy - forcefieldEnergy,
								conf.getScore(),
								breakdown.epmol.energy,
								conf.getScore() - breakdown.epmol.energy
							));
						}
					}

					for (SimpleConfSpace.Position pos2 : confSpace.positions) {

						if (pos2.index >= pos1.index) {
							break;
						}

						double boundEnergy = scoreBreakdown.getPairwise(pos1.index, 0, pos2.index, 0);
						double forcefieldEnergy = forcefieldBreakdown.getPairwise(pos1.index, 0, pos2.index, 0);
						if (boundEnergy > forcefieldEnergy + errorThreshold) {
							reports.add(String.format("%-30s    cell-score %12.6f    Sequence: %s    cell-energy %12.6f    gap: %12.6f    conf-score: %12.6f    conf-energy: %12.6f    gap: %12.6f",
								"Pair: " + pos1.formatConfPos(conf) + ":" + pos2.formatConfPos(conf),
								boundEnergy,
								confSpace.formatConf(conf),
								forcefieldEnergy,
								boundEnergy - forcefieldEnergy,
								conf.getScore(),
								breakdown.epmol.energy,
								conf.getScore() - breakdown.epmol.energy
							));
						}
					}
				}

				return reports;

			}, (reports) -> {

				if (!reports.isEmpty()) {
					try (FileWriter out = new FileWriter(outFile, true)) {
						for (String report : reports) {
							out.write(report);
							out.write("\n");
						}
					} catch (IOException ex) {
						throw new Error(ex);
					}
				}

				progress.incrementProgress();
			});
		}

		confEcalc.tasks.waitForFinish();
	}
}
