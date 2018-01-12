package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.TimeTools;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;


public class NewMethPlayground {

	public static final MathContext decimalPrecision = new MathContext(64, RoundingMode.HALF_UP);
	public static final BoltzmannCalculator bcalc = new BoltzmannCalculator(decimalPrecision);

	public static void main(String[] args) {

		// TODO: more sequences!
		// make a simple conf space with two sequences
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "THR", "PHE", "ARG").addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A6").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "THR", "PHE", "ARG").addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A7").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "THR", "PHE", "ARG").addWildTypeRotamers().setContinuous();

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// open the DB
		File confDBFile = new File("newMeth.conf.db");
		if (confDBFile.exists()) {
			confDBFile.delete();
		}
		ConfDB db = new ConfDB(confSpace, confDBFile);

		new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.use((ecalc) -> {

				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
					.build()
					.calcEnergyMatrix();

				ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
					.setTraditional()
					.setShowProgress(true)
					.build();

				log("total A* confs: %s", formatBig(astar.getNumConformations()));

				// count the number of sequences
				BigInteger numSequences = BigInteger.ZERO;
				for (SimpleConfSpace.Position pos : confSpace.positions) {
					BigInteger numResTypes = BigInteger.valueOf(pos.resFlex.resTypes.size());
					if (MathTools.isZero(numSequences)) {
						numSequences = numResTypes;
					} else {
						numSequences = numSequences.add(numResTypes);
					}
				}
				log("total A* sequences: %s", formatBig(numSequences));

				// compute lower bounds for the first K_1 conformations in the A* tree
				long maxNumConfsLowerBounded = 1000L;
				long numConfsLowerBounded = 0L;
				while (true) {

					ConfSearch.ScoredConf conf = astar.nextConf();
					if (conf == null) {
						break;
					}

					ConfDB.SequenceDB sdb = db.getSequence(confSpace.makeSequenceFromConf(conf));

					sdb.setLowerBound(
						conf.getAssignments(),
						conf.getScore(),
						TimeTools.getTimestampNs()
					);
					sdb.updateLowerEnergyOfUnsampledConfs(conf.getScore());

					if (++numConfsLowerBounded >= maxNumConfsLowerBounded) {
						break;
					}
				}

				log("Finished lower bounding conf energies");

				// sort sequences by pfunc descending upper bounds
				List<SequenceInfo> sequencesByPfuncUpperBound = new ArrayList<>();
				for (Sequence sequence : db.getSequences()) {
					SequenceInfo info = new SequenceInfo(sequence);
					info.estimatePfuncUpperBound(db.getSequence(sequence));
					sequencesByPfuncUpperBound.add(info);
				}
				sequencesByPfuncUpperBound.sort((a, b) -> b.pfuncUpperBound.compareTo(a.pfuncUpperBound));

				// pick the best K_2 sequences by pfunc upper bounds
				int maxNumBestSequences = 5;
				int numBestSequences = Math.min(maxNumBestSequences, sequencesByPfuncUpperBound.size());
				List<SequenceInfo> bestSequences = sequencesByPfuncUpperBound.subList(0, numBestSequences);
				log("%d/%s best sequences by PFunc upper bounds:", maxNumBestSequences, formatBig(numSequences));
				for (SequenceInfo info : bestSequences) {
					log("Best Sequences: %s -> %e", info.sequence, info.pfuncUpperBound);
				}

				// estimate lower bounds for the best sequences using at most K_3 confs
				int maxNumConfsUpperBounded = 100;
				for (SequenceInfo info : bestSequences) {
					log("Estimating PFunc lower bound for %s ...", info.sequence);
					info.estimatePfuncLowerBound(confEcalc, db.getSequence(info.sequence), maxNumConfsUpperBounded);
					System.out.println("\n\n" + info.toString());
				}

			}); // ecalc

		// TODO: separate energy calculation from PFunc stuff
		// TODO: visualize sequence PFunc values
	}

	public static class SequenceInfo {

		public final Sequence sequence;

		public BigInteger numConfs;
		public long numConfsLowerBounded;
		public BigInteger numConfsNotLowerBounded;
		public long numConfsUpperBounded;

		public BigDecimal pfuncUpperBoundSampled = null;
		public BigDecimal pfuncUpperBoundUnsampled = null;
		public BigDecimal pfuncUpperBound = null;
		public BigDecimal pfuncLowerBound = null;

		public SequenceInfo(Sequence sequence) {

			this.sequence = sequence;

			// count the number of confs
			numConfs = BigInteger.ZERO;
			for (SimpleConfSpace.Position pos : sequence.confSpace.positions) {
				String resType = sequence.get(pos);
				BigInteger numResConfs = BigInteger.valueOf(pos.resConfs.stream()
					.filter((resConf) -> resConf.template.name.equals(resType))
					.count());

				if (MathTools.isZero(numConfs)) {
					numConfs = numResConfs;
				} else {
					numConfs = numConfs.multiply(numResConfs);
				}
			}

			numConfsLowerBounded = 0;
			numConfsNotLowerBounded = numConfs;
			numConfsUpperBounded = 0;
		}

		public void estimatePfuncUpperBound(ConfDB.SequenceDB sdb) {

			pfuncUpperBoundSampled = BigDecimal.ZERO;
			numConfsLowerBounded = 0;
			for (ConfDB.Conf conf : sdb) {
				pfuncUpperBoundSampled = MathTools.bigAdd(
					pfuncUpperBoundSampled,
					bcalc.calc(conf.lower.energy),
					decimalPrecision
				);
				numConfsLowerBounded++;
			}
			numConfsNotLowerBounded = numConfs.subtract(BigInteger.valueOf(numConfsLowerBounded));
			pfuncUpperBoundUnsampled = bcalc.calc(sdb.getLowerEnergyOfUnsampledConfs()).multiply(new BigDecimal(numConfsNotLowerBounded), decimalPrecision);
			pfuncUpperBound = pfuncUpperBoundSampled.add(pfuncUpperBoundUnsampled, decimalPrecision);
		}

		@Override
		public String toString() {
			StringBuilder buf = new StringBuilder();
			log(buf, "sequence: %s", sequence);
			log(buf, "\ttotal confs:                   %s", formatBig(numConfs));
			log(buf, "\tnum confs lower bounded:       %d", numConfsLowerBounded);
			log(buf, "\tnum confs not lower bounded:   %s", formatBig(numConfsNotLowerBounded));
			log(buf, "\tnum confs upper bounded:       %d", numConfsUpperBounded);
			if (pfuncUpperBound != null) {
				log(buf, "\tpfunc upper bound:             %s", formatBig(pfuncUpperBound));
				log(buf, "\t\tfrom sampled:                    %s", formatBig(pfuncUpperBoundSampled));
				log(buf, "\t\tfrom unsampled:                  %s", formatBig(pfuncUpperBoundUnsampled));
				log(buf, "\t\tpercent sampled:                 %.4f", MathTools.bigDivide(pfuncUpperBoundSampled, pfuncUpperBound, decimalPrecision).doubleValue()*100.0f);
			}
			if (pfuncLowerBound != null) {
				log(buf, "\tpfunc lower bound:             %s", formatBig(pfuncLowerBound));
				double logUncertainty = Math.log10(pfuncUpperBound.doubleValue()) - Math.log10(pfuncLowerBound.doubleValue());
				log(buf, "\tpfunc log uncertainty:         %.2f", logUncertainty);
			}
			return buf.toString();
		}

		public void estimatePfuncLowerBound(ConfEnergyCalculator confEcalc, ConfDB.SequenceDB sdb, int maxNumConfs) {

			// sort confs by lower bound
			List<ConfDB.Conf> confs = new ArrayList<>();
			for (ConfDB.Conf conf : sdb) {
				confs.add(conf);
			}
			confs.sort(Comparator.comparingDouble((conf) -> conf.lower.energy));

			// TODO: use time-based limits instead of conf-based limits?
			// TODO: use smart early exit of some kind?

			// estimate the pfunc lower bound
			pfuncLowerBound = BigDecimal.ZERO;
			maxNumConfs = Math.min(maxNumConfs, confs.size());
			for (ConfDB.Conf conf : confs.subList(0, maxNumConfs)) {
				RCTuple frag = new RCTuple(conf.assignments);
				confEcalc.calcEnergyAsync(frag, (epmol) -> {

					// NOTE: this is on a listener thread, separate from the main thread

					sdb.setUpperBound(
						conf.assignments,
						epmol.energy,
						TimeTools.getTimestampNs()
					);

					pfuncLowerBound = MathTools.bigAdd(
						pfuncLowerBound,
						bcalc.calc(epmol.energy),
						decimalPrecision
					);

					numConfsUpperBounded++;
				});
			}

			confEcalc.ecalc.tasks.waitForFinish();
		}
	}

	private static void log(String format, Object ... args) {
		System.out.println(String.format(format, args));
	}

	private static void log(StringBuilder buf, String format, Object ... args) {
		buf.append(String.format(format + "\n", args));
	}

	private static String formatBig(BigInteger i) {
		if (i.compareTo(BigInteger.valueOf(1000000)) < 0) {
			return String.format("%s", i);
		} else {
			return String.format("%e", i.doubleValue());
		}
	}

	private static String formatBig(BigDecimal f) {
		return String.format("%e (%.2f)", f.doubleValue(), Math.log10(f.doubleValue()));
	}
}
