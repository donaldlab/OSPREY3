package edu.duke.cs.osprey.kstar.pfunc;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayDeque;
import java.util.Deque;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.InvertedPruningMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;

public class SimplePartitionFunction implements PartitionFunction {
	
    private EnergyMatrix emat;
    private PruningMatrix pmat;
    private ConfSearchFactory confSearchFactory;
    private ConfEnergyCalculator.Async ecalc;
    
	private double targetEpsilon;
	private Status status;
	private Values values;
	private ConfSearch tree;
	private BoltzmannCalculator boltzmann;
	private Deque<ScoredConf> confs;
	private int numConfsEvaluated;
	private BigInteger numConfsRemaining;
	private BigDecimal sumScoreWeights;
	private BigDecimal boundOnRemainingScoreWeights;
	private Stopwatch stopwatch;
	
	public SimplePartitionFunction(EnergyMatrix emat, PruningMatrix pmat, ConfSearchFactory confSearchFactory, ConfEnergyCalculator.Async ecalc) {
		
		this.emat = emat;
		this.pmat = pmat;
		this.confSearchFactory = confSearchFactory;
		this.ecalc = ecalc;
		
		targetEpsilon = Double.NaN;
		status = null;
		values = null;
		tree = null;
		boltzmann = new BoltzmannCalculator();
		confs = null;
		numConfsEvaluated = 0;
		numConfsRemaining = null;
		sumScoreWeights = null;
		boundOnRemainingScoreWeights = null;
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
	public void init(double targetEpsilon) {
		
		this.targetEpsilon = targetEpsilon;
		
		status = Status.Estimating;
		values = new Values();
		
		// compute p*: boltzmann-weight the scores for all pruned conformations
		values.pstar = calcWeightSumUpperBound(confSearchFactory.make(emat, new InvertedPruningMatrix(pmat)));
		
		// make the search tree for computing q*
		tree = confSearchFactory.make(emat, pmat);
		confs = new ArrayDeque<>();
		numConfsEvaluated = 0;
		numConfsRemaining = tree.getNumConformations();
		sumScoreWeights = BigDecimal.ZERO;
		boundOnRemainingScoreWeights = BigDecimal.ZERO;
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
			
			// report progress
			// TODO: make configurable?
			System.out.println(String.format("remaining: %d, sum: %10f, bound-remaining: %10f, bound-all: %10f, epsilon: %10f",
				numConfsRemaining, sum, boundOnRemaining, boundOnAll, effectiveEpsilon
			));
			
			if (effectiveEpsilon <= targetEpsilon) {
				break;
			}
		}
		
		return boundOnAll;
	}

	@Override
	public void step() {
		
		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}
		
		// precompute bounds on q' first by looking ahead in the conf tree
		while (true) {
			
			// read a conf from the tree
			ScoredConf conf = tree.nextConf();
			if (conf == null) {
				break;
			}
			
			numConfsRemaining = numConfsRemaining.subtract(BigInteger.ONE);
			
			BigDecimal scoreWeight = boltzmann.calc(conf.getScore());
			if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
				break;
			}
			
			sumScoreWeights = sumScoreWeights.add(scoreWeight);
			boundOnRemainingScoreWeights = scoreWeight.multiply(new BigDecimal(numConfsRemaining));
			
			// save the confs for later so we can update q*
			confs.add(conf);
			
			// stop if the bound is tight enough
			BigDecimal boundOnAll = sumScoreWeights.add(boundOnRemainingScoreWeights);
			double effectiveEpsilon = boundOnRemainingScoreWeights.divide(boundOnAll, RoundingMode.HALF_UP).doubleValue();
			if (effectiveEpsilon <= 0.01) {
				break;
			}
		}
		
		// then update q*
		
		// get a conf from the queue
		if (confs.isEmpty()) {
			status = Status.NotEnoughConformations;
			return;
		}
		ScoredConf conf = confs.removeFirst();
		
		// TODO: parallel calculation
		
		// get its energy
		EnergiedConf econf = ecalc.calcEnergy(conf);
		numConfsEvaluated++;
		
		// get the boltzmann weight
		BigDecimal weight = boltzmann.calc(econf.getEnergy());
		if (weight.compareTo(BigDecimal.ZERO) == 0) {
			status = Status.NotEnoughFiniteEnergies;
			return;
		}
		
		// update q*
		values.qstar = values.qstar.add(weight);
		
		// update q'
		sumScoreWeights = sumScoreWeights.subtract(boltzmann.calc(conf.getScore()));
		values.qprime = sumScoreWeights.add(boundOnRemainingScoreWeights);
		
		// report progress
		MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
		System.out.println(String.format("conf: %4d, energy: %.6f, q*: %12e, q': %12e, epsilon: %.6f, time: %10s, heapMem: %.0f%%",
			numConfsEvaluated, econf.getEnergy(), values.qstar, values.qprime, values.getEffectiveEpsilon(),
			stopwatch.getTime(2),
			100f*heapMem.getUsed()/heapMem.getMax()
		));
		
		if (values.getEffectiveEpsilon() <= targetEpsilon) {
			status = Status.Estimated;
		}
	}
}
