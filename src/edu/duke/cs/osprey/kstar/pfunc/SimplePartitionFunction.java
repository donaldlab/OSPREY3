package edu.duke.cs.osprey.kstar.pfunc;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.pruning.InvertedPruningMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;

// NOTE: this class is mostly a copy of ParallelConfPartitionFunction
// sadly, we can't change ParallelConfPartitionFunction because we need to maintain backwards compatibility
// so the changes needed to enable Python scripting will go here instead
public class SimplePartitionFunction implements PartitionFunction {

    public final EnergyMatrix emat;
    public final PruningMatrix pmat;
    public final ConfSearchFactory confSearchFactory;
    public final ConfEnergyCalculator ecalc;

	private double targetEpsilon = Double.NaN;
	private Status status = null;
	private Values values = null;
	private BoltzmannCalculator boltzmann = new BoltzmannCalculator();
	private ConfSearch.Splitter.Stream scoreConfs = null;
	private ConfSearch.Splitter.Stream energyConfs = null;
	private int numConfsEvaluated = 0;
	private BigInteger numConfsToScore = null;
	private BigDecimal qprimeUnevaluated = null;
	private BigDecimal qprimeUnscored = null;
	private Stopwatch stopwatch = null;
	private boolean isReportingProgress = true;
	private ConfListener confListener = null;

	public SimplePartitionFunction(EnergyMatrix emat, PruningMatrix pmat, ConfSearchFactory confSearchFactory, ConfEnergyCalculator ecalc) {
		this.emat = emat;
		this.pmat = pmat;
		this.confSearchFactory = confSearchFactory;
		this.ecalc = ecalc;
	}
	
	@Override
	public void setReportProgress(boolean val) {
		isReportingProgress = val;
	}
	
	@Override
	public void setConfListener(ConfListener val) {
		confListener = val;
	}
	
	@Override
	public Status getStatus() {
		return status;
	}
	
	@Override
	public Values getValues() {
		return values;
	}
	
	@Override
	public int getNumConfsEvaluated() {
		return numConfsEvaluated;
	}
	
	@Override
	public int getParallelism() {
		return ecalc.tasks.getParallelism();
	}

	@Override
	public void init(double targetEpsilon) {
		
		this.targetEpsilon = targetEpsilon;
		
		status = Status.Estimating;
		values = new Values();
		
		// compute p*: boltzmann-weight the scores for all pruned conformations
		values.pstar = calcWeightSumUpperBound(confSearchFactory.make(emat, new InvertedPruningMatrix(pmat)));
		
		// make the search tree for computing q*
		ConfSearch tree = confSearchFactory.make(emat, pmat);

		ConfSearch.Splitter confsSplitter = new ConfSearch.Splitter(tree);
		scoreConfs = confsSplitter.makeStream();
		energyConfs = confsSplitter.makeStream();
		numConfsEvaluated = 0;
		numConfsToScore = tree.getNumConformations();
		qprimeUnevaluated = BigDecimal.ZERO;
		qprimeUnscored = BigDecimal.ZERO;
		stopwatch = new Stopwatch().start();
	}

	private BigDecimal calcWeightSumUpperBound(ConfSearch tree) {
		
		BigDecimal sum = BigDecimal.ZERO;
		BigDecimal boundOnAll = BigDecimal.ZERO;
		
		BigInteger numConfsRemaining = tree.getNumConformations();
		
		while (true) {
			
			// get the next conf
			ScoredConf conf = tree.nextConf();
			if (conf == null) {
				break;
			}
			
			// compute the boltzmann weight for this conf
			BigDecimal weight = boltzmann.calc(conf.getScore());
			if (weight.compareTo(BigDecimal.ZERO) == 0) {
				break;
			}
			
			// update the sum
			sum = sum.add(weight);
			
			// update the upper bound on the remaining sum
			numConfsRemaining = numConfsRemaining.subtract(BigInteger.ONE);
			BigDecimal boundOnRemaining = weight.multiply(new BigDecimal(numConfsRemaining));
			
			// update the upper bound on the total sum
			boundOnAll = sum.add(boundOnRemaining);
			
			// stop if the bound is tight enough
			double effectiveEpsilon = boundOnRemaining.divide(boundOnAll, RoundingMode.HALF_UP).doubleValue();
			if (effectiveEpsilon <= 0.01) {
				break;
			}
		}
		
		return boundOnAll;
	}
	
	@Override
	public void compute() {
		compute(Integer.MAX_VALUE);
	}

	@Override
	public void compute(int maxNumConfs) {

		if (status == null) {
			throw new IllegalStateException("pfunc was not initialized. Call init() before compute()");
		}
		
		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}
		
		int stopAtConf = numConfsEvaluated + maxNumConfs;
		while (true) {
			
			ScoredConf conf;
			
			// sync to keep from racing the listener thread
			synchronized (this) {
			
				// did we win?
				boolean hitEpsilonTarget = values.getEffectiveEpsilon() <= targetEpsilon;
				if (hitEpsilonTarget) {
					status = Status.Estimated;
					break;
				}
				
				// should we stop anyway?
				if (!status.canContinue() || numConfsEvaluated >= stopAtConf) {
					break;
				}
				
				// try another conf
				conf = energyConfs.nextConf();
				if (conf == null) {
					status = Status.OutOfConformations;
				} else if (Double.isInfinite(conf.getScore()) || boltzmann.calc(conf.getScore()).compareTo(BigDecimal.ZERO) == 0) {
					status = Status.OutOfLowEnergies;
				} else {
					status = Status.Estimating;
				}
			}
			
			if (!status.canContinue()) {
				// we hit a failure condition and need to stop
				// but wait for current async minimizations to finish in case that pushes over the epsilon target
				// NOTE: don't wait for finish inside the sync, since that will definitely deadlock
				ecalc.tasks.waitForFinish();
				continue;
			}
			
			// do the energy calculation asynchronously
			ecalc.calcEnergyAsync(conf, (EnergiedConf econf) -> {

				// we're on the listener thread, so sync to keep from racing the main thread
				synchronized (SimplePartitionFunction.this) {
					
					// energy calculation done
					// if this conf has infinite energy, just ignore it
					// another later conf could still have finite energy
					if (Double.isInfinite(econf.getEnergy())) {
						return;
					}
					
					// update pfunc state
					values.qstar = values.qstar.add(boltzmann.calc(econf.getEnergy()));
					values.qprime = updateQprime(econf);
					
					// report progress if needed
					if (isReportingProgress) {
						MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
						System.out.println(String.format("conf: %4d, score: %.6f, energy: %.6f, q*: %12e, q': %12e, epsilon: %.6f, time: %10s, heapMem: %.0f%%",
							numConfsEvaluated, econf.getScore(), econf.getEnergy(), values.qstar, values.qprime, values.getEffectiveEpsilon(),
							stopwatch.getTime(2),
							100f*heapMem.getUsed()/heapMem.getMax()
						));
					}
					
					// report confs if needed
					if (confListener != null) {
						confListener.onConf(econf);
					}
					
					numConfsEvaluated++;
				}
			});
		}
		
		// wait for any remaining async minimizations to finish
		ecalc.tasks.waitForFinish();
	}

	private BigDecimal updateQprime(EnergiedConf econf) {
		
		// look through the conf tree to get conf scores
		// (which should be lower bounds on the conf energy)
		while (true) {
			
			// read a conf from the tree
			ScoredConf conf = scoreConfs.nextConf();
			if (conf == null) {
				break;
			}
			
			// get the boltzmann weight
			BigDecimal scoreWeight = boltzmann.calc(conf.getScore());
			if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
				break;
			}
			
			// update q' parts
			numConfsToScore = numConfsToScore.subtract(BigInteger.ONE);
			qprimeUnevaluated = qprimeUnevaluated.add(scoreWeight);
			qprimeUnscored = scoreWeight.multiply(new BigDecimal(numConfsToScore));
			
			// stop if the bound on q' is tight enough
			double tightness = qprimeUnscored.divide(qprimeUnevaluated.add(qprimeUnscored), RoundingMode.HALF_UP).doubleValue();
			if (tightness <= 0.01) {
				break;
			}
		}
		
		qprimeUnevaluated = qprimeUnevaluated.subtract(boltzmann.calc(econf.getScore()));
		return qprimeUnevaluated.add(qprimeUnscored);
	}
}
