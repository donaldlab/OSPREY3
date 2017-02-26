package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.ParallelConfPartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;

public class ParallelPartitionFunction extends ParallelConfPartitionFunction {

	protected BigDecimal qstarScoreWeights;
	int numActiveThreads;

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
		numActiveThreads = 0;
		stopwatch = new Stopwatch().start();
	}
	
	@Override
	protected BigDecimal updateQprime(EnergiedConf econf) {
		
		// look through the conf tree to get conf scores
		// (which should be lower bounds on the conf energy)
		while (true) {
			
			// read a conf from the tree
			ScoredConf conf = scoreConfs.next();
			if (conf == null) {
				qprimeUnscored = BigDecimal.ZERO;
				break;
			}
			
			// get the boltzmann weight
			BigDecimal scoreWeight = boltzmann.calc(conf.getScore());
			if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
				qprimeUnscored = BigDecimal.ZERO;
				break;
			}
			
			// update q' parts
			numConfsToScore = numConfsToScore.subtract(BigInteger.ONE);
			qprimeUnevaluated = qprimeUnevaluated.add(scoreWeight);
			qprimeUnscored = scoreWeight.multiply(new BigDecimal(numConfsToScore));
			
			// stop if the bound on q' is tight enough
			double effectiveEpsilon = qprimeUnscored.divide(qprimeUnevaluated.add(qprimeUnscored), RoundingMode.HALF_UP).doubleValue();
			if (effectiveEpsilon <= 0.01) {
				break;
			}
		}
		
		qprimeUnevaluated = qprimeUnevaluated.subtract(boltzmann.calc(econf.getScore()));
		return qprimeUnevaluated.add(qprimeUnscored);
	}

	@Override
	public void compute(int maxNumConfs) {
		numActiveThreads = 0;
		
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

				if ((conf = energyConfs.next()) == null) {
					while(numActiveThreads > 0) {
						try { this.wait(); } catch (InterruptedException e) { e.printStackTrace(); }
					}
					if(status != Status.Estimated) status = Status.NotEnoughConformations;
					break;
				}

				if (boltzmann.calc(conf.getScore()).compareTo(BigDecimal.ZERO) == 0) {
					while(numActiveThreads > 0) {
						try { this.wait(); } catch (InterruptedException e) { e.printStackTrace(); }
					}
					if(status != Status.Estimated) status = Status.NotEnoughFiniteEnergies;
					break;
				}
				
				++numActiveThreads;
			}

			// do the energy calculation asynchronously
			ecalc.calcEnergyAsync(conf, (EnergiedConf econf) -> {

				// energy calculation done

				// this is (potentially) running on a task executor listener thread
				// so lock to keep from racing the main thread
				synchronized (ParallelPartitionFunction.this) {

					// get the boltzmann weight
					BigDecimal energyWeight = boltzmann.calc(econf.getEnergy());

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

					--numActiveThreads;
					this.notify();
				}
			});
		}
		
		// wait for any remaining async minimizations to finish
		ecalc.waitForFinish();
	}

	public void compute(BigDecimal targetScoreWeights) {
		numActiveThreads = 0;
		
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
				if (!status.canContinue() || qstarScoreWeights.compareTo(targetScoreWeights) >= 0) {
					break;
				}

				if ((conf = energyConfs.next()) == null) {
					while(numActiveThreads > 0) {
						try { this.wait(); } catch (InterruptedException e) { e.printStackTrace(); }
					}
					if(status != Status.Estimated) status = Status.NotEnoughConformations;
					break;
				}

				if (boltzmann.calc(conf.getScore()).compareTo(BigDecimal.ZERO) == 0) {
					while(numActiveThreads > 0) {
						try { this.wait(); } catch (InterruptedException e) { e.printStackTrace(); }
					}
					if(status != Status.Estimated) status = Status.NotEnoughFiniteEnergies;
					break;
				}
				
				++numActiveThreads;
			}

			// do the energy calculation asynchronously
			ecalc.calcEnergyAsync(conf, (EnergiedConf econf) -> {

				// energy calculation done

				// this is (potentially) running on a task executor listener thread
				// so lock to keep from racing the main thread
				synchronized (ParallelPartitionFunction.this) {

					// get the boltzmann weight
					BigDecimal scoreWeight = boltzmann.calc(econf.getScore());
					qstarScoreWeights = qstarScoreWeights.add(scoreWeight);	
					BigDecimal energyWeight = boltzmann.calc(econf.getEnergy());

					// update pfunc state
					numConfsEvaluated++;
					values.qstar = values.qstar.add(energyWeight);
					values.qprime = updateQprime(econf);
					BigDecimal pdiff = targetScoreWeights.subtract(qstarScoreWeights);

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
					
					--numActiveThreads;
					this.notify();
				}
			});
		}

		// wait for any remaining async minimizations to finish
		ecalc.waitForFinish();
	}

}
