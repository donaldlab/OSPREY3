package edu.duke.cs.osprey.kstar.pfunc;

import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.math.BigDecimal;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.multistatekstar.QPruningMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;

public class ParallelPartitionFunction extends ParallelConfPartitionFunction {

	protected BigDecimal qstarScoreWeights;
	
	public ParallelPartitionFunction(EnergyMatrix emat, PruningMatrix pmat, ConfSearchFactory confSearchFactory,
			Async ecalc) {
		super(emat, pmat, confSearchFactory, ecalc);
	}

	@Override
	public void init(double targetEpsilon) {

		this.targetEpsilon = targetEpsilon;

		status = Status.Estimating;
		values = new Values();

		// compute p*: boltzmann-weight the scores for all pruned conformations
		values.pstar = calcWeightSumUpperBound(confSearchFactory.make(emat, ((QPruningMatrix)pmat).invert()));

		// make the search tree for computing q*
		ConfSearch tree = confSearchFactory.make(emat, pmat);
		ConfSearch.Splitter confsSplitter = new ConfSearch.Splitter(tree);
		scoreConfs = confsSplitter.makeStream();
		energyConfs = confsSplitter.makeStream();
		numConfsEvaluated = 0;
		numConfsToScore = tree.getNumConformations();
		qprimeUnevaluated = BigDecimal.ZERO;
		qprimeUnscored = BigDecimal.ZERO;
		qstarScoreWeights = BigDecimal.ZERO;
		stopwatch = new Stopwatch().start();
	}
	
	@Override
	public void compute(int maxNumConfs) {
		
		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}
		
		int stopAtConf = numConfsEvaluated + maxNumConfs;
		while (true) {
			
			// wait for space to open up Before getting a new conf to minimize
			ecalc.waitForSpace();
			
			// get a conf from the tree
			// lock though to keep from racing the listener thread on the conf tree
			ScoredConf conf;
			synchronized (this) {
				
				// should we keep going?
				if (!status.canContinue() || numConfsEvaluated >= stopAtConf) {
					break;
				}
				
				conf = energyConfs.next();
				if (conf == null) {
					status = Status.NotEnoughConformations;
					return;
				}
			}
			
			// do the energy calculation asynchronously
			ecalc.calcEnergyAsync(conf, (EnergiedConf econf) -> {
					
				// energy calculation done
				
				// this is (potentially) running on a task executor listener thread
				// so lock to keep from racing the main thread
				synchronized (ParallelPartitionFunction.this) {
				
					// get the boltzmann weight
					BigDecimal energyWeight = boltzmann.calc(econf.getEnergy());
					BigDecimal scoreWeight = boltzmann.calc(econf.getScore());
					if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
						status = Status.NotEnoughFiniteEnergies;
						return;
					}
					
					// update pfunc state
					numConfsEvaluated++;
					values.qstar = values.qstar.add(energyWeight);
					values.qprime = updateQprime(econf);
					
					// report progress if needed
					if (isReportingProgress) {
						MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
						System.out.println(String.format("conf: %4d, energy: %.6f, q*: %12e, q': %12e, p*: %12e, epsilon: %.6f, time: %10s, heapMem: %.0f%%",
							numConfsEvaluated, econf.getEnergy(), values.qstar, values.qprime, values.pstar, values.getEffectiveEpsilon(),
							stopwatch.getTime(2),
							100f*heapMem.getUsed()/heapMem.getMax()
						));
					}
					
					// report confs if needed
					if (confListener != null) {
						confListener.onConf(econf);
					}
					
					// update status if needed
					if (values.getEffectiveEpsilon() <= targetEpsilon) {
						status = Status.Estimated;
					}
				}
			});
		}
		
		// wait for any remaining async minimizations to finish
		ecalc.waitForFinish();
	}
	
	public void compute(BigDecimal target) {

		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		while (true) {

			// wait for space to open up Before getting a new conf to minimize
			ecalc.waitForSpace();

			// get a conf from the tree
			// lock though to keep from racing the listener thread on the conf tree
			ScoredConf conf;
			synchronized (this) {

				// should we keep going?
				if (!status.canContinue() || qstarScoreWeights.compareTo(target) >= 0) {
					break;
				}

				conf = energyConfs.next();
				if (conf == null) {
					status = Status.NotEnoughConformations;
					break;
				}
			}

			// do the energy calculation asynchronously
			ecalc.calcEnergyAsync(conf, (EnergiedConf econf) -> {

				// energy calculation done

				// this is (potentially) running on a task executor listener thread
				// so lock to keep from racing the main thread
				synchronized (ParallelPartitionFunction.this) {

					// get the boltzmann weight
					BigDecimal energyWeight = boltzmann.calc(econf.getEnergy());
					BigDecimal scoreWeight = boltzmann.calc(econf.getScore());
					if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
						status = Status.NotEnoughFiniteEnergies;
						return;
					}
					
					qstarScoreWeights = qstarScoreWeights.add(scoreWeight);

					// update pfunc state
					numConfsEvaluated++;
					values.qstar = values.qstar.add(energyWeight);
					values.qprime = updateQprime(econf);
					BigDecimal pdiff = target.subtract(qstarScoreWeights);

					// report progress if needed
					if (isReportingProgress) {
						MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
						System.out.println(String.format("conf: %4d, energy: %.6f, q*: %12e, q': %12e, p*-p1*: %12e, time: %10s, heapMem: %.0f%%",
								numConfsEvaluated, econf.getEnergy(), values.qstar, values.qprime, pdiff,
								stopwatch.getTime(2),
								100f*heapMem.getUsed()/heapMem.getMax()
								));
					}

					// report confs if needed
					if (confListener != null) {
						confListener.onConf(econf);
					}

					// update status if needed
					if (values.getEffectiveEpsilon() <= targetEpsilon) {
						status = Status.Estimated;
					}
				}
			});
		}

		// wait for any remaining async minimizations to finish
		ecalc.waitForFinish();
	}

}
