package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfRanker;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;

import java.io.File;
import java.math.BigInteger;
import java.util.function.Consumer;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;


public class BenchmarkConfRanker {

	public static void main(String[] args) {

		// try that gp120 design with massive flexibility
		Molecule mol = PDBIO.readResource("/gp120/gp120SRVRC26.09SR.pdb");

		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder()
			.addMoleculeForWildTypeRotamers(mol)
			.addTemplates(FileTools.readResource("/gp120/all_nuc94_and_gr.in"))
			.addTemplateCoords(FileTools.readResource("/gp120/all_amino_coords.in"))
			.addRotamers(FileTools.readResource("/gp120/GenericRotamers.dat"))
			.build();

		Strand ligand = new Strand.Builder(mol)
			.setResidues("H1792", "L2250")
			.setTemplateLibrary(templateLib)
			.build();
		ligand.flexibility.get("H1901").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1904").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIS", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1905").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIS", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1906").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1907").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1908").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		Strand target = new Strand.Builder(mol)
			.setResidues("F379", "J1791")
			.setTemplateLibrary(templateLib)
			.build();
		target.flexibility.get("G973").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G977").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G978").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G979").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G980").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("J1448").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.addStrand(target)
			.build();

		// calc the emat
		EnergyMatrix emat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()
		) {
			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
			emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File("emat.gp120.complex.dat"))
				.build()
				.calcEnergyMatrix();
		}
		// pick the wild-type sequence
		Sequence sequence = confSpace.makeWildTypeSequence();
		RCs rcs = sequence.makeRCs();

		log("confs: %s", formatBig(makeAStar(emat, rcs).getNumConformations()));

		// get the min,max scores
		double minScore = makeAStar(emat, rcs)
			.nextConf()
			.getScore();
		log("min score: %.4f", minScore);
		double maxScore = -makeAStar(new NegatedEnergyMatrix(confSpace, emat), rcs)
			.nextConf()
			.getScore();
		log("max score: %.4f", maxScore);

		// pick a few energy thresholds to rank, relative to the min score
		double[] energyOffsets = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };//, 11, 12 };//, 13, 14 };

		log("");

		// how quickly can we just enumerate the confs the regular way?
		{
			int offsetIndex = 0;
			long numConfs = 0;
			long[] times = new long[energyOffsets.length];
			ConfAStarTree astar = makeAStar(emat, rcs);
			Stopwatch stopwatch = new Stopwatch().start();
			while (offsetIndex < energyOffsets.length) {

				ConfSearch.ScoredConf conf = astar.nextConf();
				if (conf == null) {
					break;
				}

				numConfs++;

				if (conf.getScore() >= minScore + energyOffsets[offsetIndex]) {
					times[offsetIndex] = stopwatch.getTimeNs();
					log("energy %.4f has %12d confs in %10s  (%10s more than last)",
						conf.getScore(), numConfs, TimeFormatter.format(times[offsetIndex], 2),
						offsetIndex > 0 ? TimeFormatter.format(times[offsetIndex] - times[offsetIndex - 1], 2) : "no"
					);
					offsetIndex++;
				}
			}
		}

		// how quickly can we rank the confs?
		ConfRanker ranker = new ConfRanker.Builder(confSpace, emat)
			.setRCs(rcs)
			//.setReportProgress(true)
			.build();

		Consumer<Double> func = (queryScore) -> {
			Stopwatch stopwatch = new Stopwatch().start();
			BigInteger rank = ranker.getNumConfsAtMost(queryScore);
			stopwatch.stop();
			long timeNs = stopwatch.getTimeNs();
			log("energy %.4f has %12s confs in %10s  (%10s more than last)",
				queryScore, formatBig(rank),
				TimeFormatter.format(timeNs, 2),
				TimeFormatter.format(timeNs, 2)
			);
		};

		log("");

		// these should be super duper fast
		func.accept(minScore - 0.00001);
		func.accept(maxScore + 0.1);

		log("");

		for (double offset : energyOffsets) {
			func.accept(minScore + offset);
		}

		// this one takes waaaaay too long to finish
		// I let it run for like 24 hours once...
		//func.accept(0.0);
	}

	private static ConfAStarTree makeAStar(EnergyMatrix emat, RCs rcs) {
		ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
			.setTraditional()
			.build();
		astar.setParallelism(Parallelism.makeCpu(6));
		return astar;
	}
}
