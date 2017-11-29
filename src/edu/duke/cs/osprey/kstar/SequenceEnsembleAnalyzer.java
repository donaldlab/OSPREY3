package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.JvmMem;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.*;

/**
 * Shows information about a single sequence.
 */
public class SequenceEnsembleAnalyzer {

	public class Analysis {

		public final KStar.ConfSpaceType type;
		public final KStar.Sequence complexSequence;
		public final KStar.Sequence filteredSequence;
		public final ConfAnalyzer.EnsembleAnalysis ensemble;

		public Analysis(KStar.ConfSpaceType type, KStar.Sequence complexSequence, KStar.Sequence filteredSequence, ConfAnalyzer.EnsembleAnalysis ensemble) {
			this.type = type;
			this.complexSequence = complexSequence;
			this.filteredSequence = filteredSequence;
			this.ensemble = ensemble;
		}

		public SimpleConfSpace getConfSpace() {
			return infosByType.get(type).confSpace;
		}

		public void writePdbs(String filePattern) {
			ensemble.writePdbs(filePattern);
		}

		@Override
		public String toString() {

			SimpleConfSpace confSpace = getConfSpace();
			int indexSize = 1 + (int)Math.log10(ensemble.analyses.size());

			StringBuilder buf = new StringBuilder();
			buf.append(String.format("Residues           %s\n", confSpace.formatResidueNumbers()));
			buf.append(String.format("%-16s   %s\n", type.name() + " Sequence", filteredSequence.toString(KStar.Sequence.makeWildType(confSpace), 5)));
			buf.append(String.format("Ensemble of %d conformations:\n", ensemble.analyses.size()));
			for (int i=0; i<ensemble.analyses.size(); i++) {
				ConfAnalyzer.ConfAnalysis analysis = ensemble.analyses.get(i);
				buf.append("\t");
				buf.append(String.format("%" + indexSize + "d/%" + indexSize + "d", i + 1, ensemble.analyses.size()));
				buf.append(String.format("     Energy: %-12.6f     Score: %-12.6f", analysis.epmol.energy, analysis.score));
				buf.append("     Rotamers: ");
				buf.append(confSpace.formatConfRotamers(analysis.assignments));
				buf.append("     Residue Conf IDs: ");
				buf.append(SimpleConfSpace.formatConfRCs(analysis.assignments));
				buf.append("\n");
			}
			return buf.toString();
		}
	}

	public class ConfSpaceInfo {

		public final KStar.ConfSpaceType type;
		public final SimpleConfSpace confSpace;
		public final ConfEnergyCalculator confEcalc;
		public final Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToThisMap;

		public EnergyMatrix emat = null;

		public ConfSpaceInfo(KStar.ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc, Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToThisMap) {
			this.type = type;
			this.confSpace = confSpace;
			this.confEcalc = confEcalc;
			this.complexToThisMap = complexToThisMap;
		}

		public void calcEmat() {
			SimplerEnergyMatrixCalculator.Builder builder = new SimplerEnergyMatrixCalculator.Builder(confEcalc);
			if (settings.energyMatrixCachePattern != null) {
				builder.setCacheFile(new File(settings.applyEnergyMatrixCachePattern(type.name().toLowerCase())));
			}
			emat = builder.build().calcEnergyMatrix();
		}

		public KStar.Sequence makeWildTypeSequence() {
			KStar.Sequence sequence = new KStar.Sequence(confSpace.positions.size());
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				sequence.set(pos.index, pos.strand.mol.getResByPDBResNumber(pos.resNum).template.name);
			}
			return sequence;
		}

		public KStar.Sequence filterSequence(KStar.Sequence complexSequence) {

			// don't need to filter complex sequences
			if (complexToThisMap == null) {
				assert (type == KStar.ConfSpaceType.Complex);
				return complexSequence;
			}

			KStar.Sequence filteredSequence = makeWildTypeSequence();
			for (SimpleConfSpace.Position pos : complex.confSpace.positions) {
				SimpleConfSpace.Position proteinPos = complexToThisMap.get(pos);
				if (proteinPos != null) {
					filteredSequence.set(proteinPos.index, complexSequence.get(pos.index));
				}
			}
			return filteredSequence;
		}
	}


	/** A configuration space containing just the protein strand */
	public final ConfSpaceInfo protein;

	/** A configuration space containing just the ligand strand */
	public final ConfSpaceInfo ligand;

	/** A configuration space containing both the protein and ligand strands */
	public final ConfSpaceInfo complex;

	/** Calculates the energy for a molecule */
	public final EnergyCalculator ecalc;

	/** A function that makes a ConfEnergyCalculator with the desired options */
	public final KStar.ConfEnergyCalculatorFactory confEcalcFactory;

	/** A function that makes a ConfSearchFactory (e.g, A* search) with the desired options */
	public final KStar.ConfSearchFactory confSearchFactory;

	/** K* settings */
	public final KStar.Settings settings;

	private final Map<KStar.ConfSpaceType,ConfSpaceInfo> infosByType;

	public SequenceEnsembleAnalyzer(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator ecalc, KStar.ConfEnergyCalculatorFactory confEcalcFactory, KStar.ConfSearchFactory confSearchFactory, KStar.Settings settings) {

		this.protein = new ConfSpaceInfo(KStar.ConfSpaceType.Protein, protein, confEcalcFactory.make(protein, ecalc), complex.mapPositionsTo(protein));
		this.ligand = new ConfSpaceInfo(KStar.ConfSpaceType.Ligand, ligand, confEcalcFactory.make(ligand, ecalc), complex.mapPositionsTo(ligand));
		this.complex = new ConfSpaceInfo(KStar.ConfSpaceType.Complex, complex, confEcalcFactory.make(complex, ecalc), null);
		this.ecalc = ecalc;
		this.confEcalcFactory = confEcalcFactory;
		this.confSearchFactory = confSearchFactory;
		this.settings = settings;

        System.out.println("Hello World!");

		// index conf space infos by type
		infosByType = new EnumMap<>(KStar.ConfSpaceType.class);
		infosByType.put(KStar.ConfSpaceType.Protein, this.protein);
		infosByType.put(KStar.ConfSpaceType.Ligand, this.ligand);
		infosByType.put(KStar.ConfSpaceType.Complex, this.complex);

		// calc emats if needed, or read them from the cache
		this.protein.calcEmat();
		this.ligand.calcEmat();
		this.complex.calcEmat();
	}

	public void analyze(KStar.Sequence complexSequence, int numConfs) {
        System.out.println("Analyzing "+complexSequence+"!");
        
	    Map<KStar.ConfSpaceType, PartitionFunction.Result> resultsByType = new EnumMap<>(KStar.ConfSpaceType.class);
		
        for(KStar.ConfSpaceType type: infosByType.keySet()) {
            resultsByType.put(type, computePFuncFor(infosByType.get(type), complexSequence, numConfs));
        }

        /* TEMP
        KStarScore score = new KStarScore(
        	resultsByType.get(KStar.ConfSpaceType.Protein),
            resultsByType.get(KStar.ConfSpaceType.Ligand),
            resultsByType.get(KStar.ConfSpaceType.Complex)
		);
        System.out.println(score.printEnsembleScoreAndBounds());
        */
        /* TODO: Write to a file, as well.
        KStar.KStarScoreWriter writer = new KStarScoreWriter.ToFile(file, new KStarScoreWriter.Formatter.Log());
        settings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
            sequenceNumber,
            n,
            complexSequence,
            complex.confSpace,
            kstarScore
        ));
        */
        
	}

    private PartitionFunction.Result computePFuncFor(ConfSpaceInfo info, KStar.Sequence complexSequence, int numConfs) {

		// filter the sequence and get the RCs
		KStar.Sequence filteredSequence = info.filterSequence(complexSequence);
		RCs rcs = filteredSequence.makeRCs(info.confSpace);

		// TEMP
		System.out.println(String.format("%s\n\t%s\n\t%s\n\tconfs: %e",
			info.type,
			complexSequence,
			filteredSequence,
			rcs.getNumConformations().doubleValue()
		));
		for (SimpleConfSpace.Position pos : info.confSpace.positions) {
			System.out.println(String.format("\t%s %s %d", pos.resNum, pos.resFlex.wildType, pos.resConfs.size()));
		}

        // make the A* search
        ConfSearch astar = confSearchFactory.make(info.emat, rcs);

		/* TEMP
        // make the partition function
        SimplePartitionFunction pfunc = new SimplePartitionFunction(astar, info.confEcalc);
        pfunc.setReportProgress(settings.showPfuncProgress);

        // compute pfunc for protein
        pfunc.init(settings.epsilon);
        pfunc.compute(numConfs);
        PartitionFunction.Result proteinResult = pfunc.makeResult();
        System.out.println(proteinResult);
        return proteinResult;
        */

		// TEMP: just compute the upper bounds
		System.out.println("computing upper bound...");
		final double upperBoundEpsilon = 0.00001;
		SimplePartitionFunction.UpperBoundCalculator upperBound = new SimplePartitionFunction.UpperBoundCalculator(astar);
		while (upperBound.delta > upperBoundEpsilon) {
			upperBound.run(10000, upperBoundEpsilon);
			System.out.println(String.format("\tconfs: %.2f M  upper bound: %e    delta: %f > %f   mem: %s",
				upperBound.numScoredConfs/1000000.0,
				upperBound.totalBound.doubleValue(),
				upperBound.delta,
				upperBoundEpsilon,
				JvmMem.getOldPool()
			));
		}
		System.out.println("uppper bound complete");

		return null;
    }
}
