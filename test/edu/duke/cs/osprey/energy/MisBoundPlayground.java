package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarNode;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.File;
import java.util.Arrays;
import java.util.List;


public class MisBoundPlayground {

	public static void main(String[] args)
	throws Exception {

		String dir = "__redacted__";

		ForcefieldParams ffparams = new ForcefieldParams();

		// load non-natural AA stuff
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.addTemplates(FileTools.readFile(dir + "/all_nuc94_and_gr.in"))
			.addTemplateCoords(FileTools.readFile(dir + "/all_amino_coords.in"))
			.addRotamers(FileTools.readFile(dir + "/GenericRotamers.dat"))
			.build();

		Strand strand = new Strand.Builder(PDBIO.readFile(dir + "/model0_antibody.pdb"))
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
			//.setShellDistance(6)
			.build();

		new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.make(1, 1, 16))
			.use((ecalc) -> {

				SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcReferenceEnergies();

				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.setReferenceEnergies(eref)
					.setEnergyPartition(EnergyPartition.Traditional)
					//.setEnergyPartition(EnergyPartition.AllOnPairs)
					.build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
					.setCacheFile(new File("emat.misbound.traditional.dat"))
					//.setCacheFile(new File("emat.misbound.allonpairs.dat"))
					.build()
					.calcEnergyMatrix();

				ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
					.setTraditional()
					.setShowProgress(true)
					.build();

				/* do the search
				new SimpleGMECFinder.Builder(astar, confEcalc)
					.build()
					.find(5.0);
				*/

				List<int[]> misboundedConfsTraditional = Arrays.asList(
					new int[] { 22, 131, 159 }, // Conformation score (-14.511145) is not a lower bound on the energy (-14.869177).
					new int[] { 27, 131, 159 } // Conformation score (-13.803699) is not a lower bound on the energy (-14.221509).
				);

				/* from Graham's script, RCs are different for some reason... maybe mutations are in different order
				List<int[]> misboundedConfsAllOnPairs = Arrays.asList(
					new int[] { 17, 135, 149 }, // Conformation score (-14.033345) is not a lower bound on the energy (-14.211877).
					new int[] { 17, 63, 149 }, // Conformation score (-13.812986) is not a lower bound on the energy (-13.971702).
					new int[] { 17, 124, 149 }, // Conformation score (-13.759669) is not a lower bound on the energy (-13.866282).
					new int[] { 17, 65, 149 }, // Conformation score (-13.691940) is not a lower bound on the energy (-13.805368).
					new int[] { 17, 126, 177 }, // Conformation score (-13.665316) is not a lower bound on the energy (-13.847284).
					new int[] { 17, 69, 149 } // Conformation score (-13.643942) is not a lower bound on the energy (-13.784523).
				);
				*/

				List<int[]> misboundedConfsAllOnPairs = Arrays.asList(
					new int[] { 22, 140, 149 }, // Conformation score (-14.033345) is not a lower bound on the energy (-14.211877).
					new int[] { 22, 68, 149 }, // Conformation score (-13.812986) is not a lower bound on the energy (-13.971702).
					new int[] { 22, 129, 149 }, // Conformation score (-13.759669) is not a lower bound on the energy (-13.866282).
					new int[] { 22, 70, 149 }, // Conformation score (-13.691940) is not a lower bound on the energy (-13.805368).
					new int[] { 22, 131, 177 }, // Conformation score (-13.665316) is not a lower bound on the energy (-13.847284).
					new int[] { 22, 74, 149 } // Conformation score (-13.643942) is not a lower bound on the energy (-13.784523).
				);

				//List<int[]> misboundedConfs = misboundedConfsTraditional;
				List<int[]> misboundedConfs = misboundedConfsAllOnPairs;

				// check the confs directly
				ConfIndex index = new ConfIndex(confSpace.positions.size());
				//for (int[] conf : misboundedConfs) {
				{ int[] conf = { 22, 131, 159 };

					// get an A* node for the conf
					LinkedConfAStarNode node = new LinkedConfAStarNode();
					for (SimpleConfSpace.Position pos : confSpace.positions) {
						node = node.assign(pos.index, conf[pos.index]);
					}

					// compute the A* g-score
					node.index(index);
					double score = astar.gscorer.calc(index, astar.rcs);

					// minimize the conf and get the energy
					double energy = confEcalc.calcEnergy(new RCTuple(conf));

					// get the sequence and rotamers
					StringBuilder sequence = new StringBuilder();
					StringBuilder rotamers = new StringBuilder();
					for (SimpleConfSpace.Position pos : confSpace.positions) {
						SimpleConfSpace.ResidueConf resConf = pos.resConfs.get(conf[pos.index]);
						if (sequence.length() > 0) {
							sequence.append(" ");
							rotamers.append(" ");
						}
						sequence.append(String.format("%3s", resConf.template.name));
						rotamers.append(String.format("%3s", pos.resConfs.get(conf[pos.index]).getRotamerCode()));
					}

					System.out.println(String.format("Mis-bounded conf:"
						+ "\n\tAssignments:   %s"
						+ "\n\tSequence:      %s"
						+ "\n\tRotamers:      %s"
						+ "\n\tScore:         %12.6f"
						+ "\n\tEnergy:        %12.6f",
						Arrays.toString(conf),
						sequence,
						rotamers,
						score,
						energy
					));

					// breakdown the forcefield
					ResidueForcefieldBreakdown.ByPosition breakdownByPos = new ResidueForcefieldBreakdown.ByPosition(confEcalc, conf);
					System.out.println("forcefield: " + breakdownByPos.breakdownForcefield(ResidueForcefieldBreakdown.Type.All));

					// breakdown the bound
					System.out.println("bound: " + breakdownByPos.breakdownBound(emat));
				}
			});
	}
}
