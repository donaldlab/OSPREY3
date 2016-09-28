package edu.duke.cs.osprey.kstar.pfunc;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.InvertedPruningMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;

public class ParallelConfPartitionFunction implements PartitionFunction {
	
    private EnergyMatrix emat;
    private PruningMatrix pmat;
    private ConfSearchFactory confSearchFactory;
    private ConfEnergyCalculator.Async ecalc;
    
	private double targetEpsilon;
	private Status status;
	private Values values;
	private BoltzmannCalculator boltzmann;
	private ConfSearch.Splitter.Stream scoreConfs;
	private ConfSearch.Splitter.Stream energyConfs;
	private int numConfsEvaluated;
	private BigInteger numConfsToScore;
	private BigDecimal qprimeUnevaluated;
	private BigDecimal qprimeUnscored;
	private Stopwatch stopwatch;
	private boolean isReportingProgress;
	
	public ParallelConfPartitionFunction(EnergyMatrix emat, PruningMatrix pmat, ConfSearchFactory confSearchFactory, ConfEnergyCalculator.Async ecalc) {
		this.emat = emat;
		this.pmat = pmat;
		this.confSearchFactory = confSearchFactory;
		this.ecalc = ecalc;
		
		targetEpsilon = Double.NaN;
		status = null;
		values = null;
		boltzmann = new BoltzmannCalculator();
		scoreConfs = null;
		energyConfs = null;
		numConfsEvaluated = 0;
		numConfsToScore = null;
		qprimeUnevaluated = null;
		qprimeUnscored = null;
		
		stopwatch = null;
		isReportingProgress = true;
	}
	
	@Override
	public void setReportProgress(boolean val) {
		isReportingProgress = val;
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
	public int getParallelism() {
		return ecalc.getParallelism();
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
			ecalc.calcEnergyAsync(conf, new ConfEnergyCalculator.Async.Listener() {
				@Override
				public void onEnergy(EnergiedConf econf) {
					
					// energy calculation done
					
					// this is (potentially) running on a task executor listener thread
					// so lock to keep from racing the main thread
					synchronized (ParallelConfPartitionFunction.this) {
					
						// do we even care about this result anymore?
						if (!status.canContinue()) {
							return;
						}
						
						// get the boltzmann weight
						BigDecimal energyWeight = boltzmann.calc(econf.getEnergy());
						if (energyWeight.compareTo(BigDecimal.ZERO) == 0) {
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
							System.out.println(String.format("conf: %4d, energy: %.6f, q*: %12e, q': %12e, epsilon: %.6f, time: %10s, heapMem: %.0f%%",
								numConfsEvaluated, econf.getEnergy(), values.qstar, values.qprime, values.getEffectiveEpsilon(),
								stopwatch.getTime(2),
								100f*heapMem.getUsed()/heapMem.getMax()
							));
						}
						
						// update status if needed
						if (values.getEffectiveEpsilon() <= targetEpsilon) {
							status = Status.Estimated;
						}
					}
				}
			});
		}
	}

	protected BigDecimal updateQprime(EnergiedConf econf) {
		
		// look through the conf tree to get conf scores
		// (which should be lower bounds on the conf energy)
		while (true) {
			
			// read a conf from the tree
			ScoredConf conf = scoreConfs.next();
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
			double effectiveEpsilon = qprimeUnscored.divide(qprimeUnevaluated.add(qprimeUnscored), RoundingMode.HALF_UP).doubleValue();
			if (effectiveEpsilon <= 0.01) {
				break;
			}
		}
		
		qprimeUnevaluated = qprimeUnevaluated.subtract(boltzmann.calc(econf.getScore()));
		return qprimeUnevaluated.add(qprimeUnscored);
	}
}
