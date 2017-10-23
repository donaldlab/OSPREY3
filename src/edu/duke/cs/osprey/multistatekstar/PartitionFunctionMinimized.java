package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;

import com.jogamp.common.util.InterruptSource.Thread;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.ParallelConfPartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class PartitionFunctionMinimized extends ParallelConfPartitionFunction {

	public static boolean SYNCHRONIZED_MINIMIZATION = false;
	public static double QPRIME_TIMEOUT_HRS = 0.0;
	protected PriorityQueue<ScoredConf> topConfs;
	protected int maxNumTopConfs;
	protected BigDecimal qstarScoreWeights;
	protected PruningMatrix invmat;
	protected ArrayList<ScoredConf> scoredConfs;
	protected ArrayList<EnergiedConf> energiedConfs;

	protected List<EnergiedConf> minGMECConfs;
	protected HashSet<ArrayList<Integer>> minGMECAssignments;
	
	protected boolean computeGMECRatio;
	protected boolean energiedGMECEnumerated;
	protected boolean scoredGMECPStarEnumerated;
	protected boolean scoredGMECQPrimeEnumerated;

	protected boolean computeMaxNumConfs;

	public PartitionFunctionMinimized(
			EnergyMatrix emat, 
			PruningMatrix pmat, 
			PruningMatrix invmat, 
			ConfSearchFactory confSearchFactory,
			Async ecalc
			) {
		super(emat, pmat, confSearchFactory, ecalc);
		this.invmat = invmat;
		this.qstarScoreWeights = null;
		this.topConfs = null;
		this.scoredConfs = null;
		this.energiedConfs = null;
		this.minGMECConfs = null;
		this.minGMECAssignments = null;
		this.computeGMECRatio = false;
		this.computeMaxNumConfs = false;
		this.boltzmann = new MSBoltzmannCalculator();
		this.maxNumTopConfs = 0;
	}

	protected void writeTopConfs(int state, MSSearchProblem search, String baseDir) {
		if(topConfs==null || topConfs.size()==0) return;
		String seq = search.settings.getFormattedSequence();
		if(isReportingProgress) {
			System.out.println("Writing top "+ topConfs.size()+" confs:");
		}
		seq = seq.replace(" ", ".");
		String dir = baseDir+File.separator+"State."+state+File.separator+seq;
		ObjectIO.makeDir(dir, false);
		for(int i=topConfs.size()-1;i>-1;--i) {
			if(isReportingProgress) {
				ScoredConf head = topConfs.peek();
				double energy = head instanceof EnergiedConf ? ((EnergiedConf)head).getEnergy() : head.getScore();
				System.out.println(String.format("conf: %4d.pdf, energy: %.6f", i, energy));
			}
			String PDBFileName = dir+File.separator+i+".pdb";
			search.outputMinimizedStruct(topConfs.poll().getAssignments(), PDBFileName);
		}
	}

	protected void saveConf(ScoredConf conf) {
		if(!topConfs.isEmpty() && topConfs.size() >= maxNumTopConfs) {

			ScoredConf head = topConfs.peek();
			double e1 = head instanceof EnergiedConf ? ((EnergiedConf)head).getEnergy() : head.getScore();
			double e2 = conf instanceof EnergiedConf ? ((EnergiedConf)conf).getEnergy() : conf.getScore();

			if(e1 > e2) topConfs.poll();
			else return;
		}
		topConfs.add(conf);
	}

	protected void saveEConfs(PriorityQueue<ScoredConf> other) {
		if(topConfs==null || other==null) return;
		while(other.size()>0) 
			saveConf(other.poll());
	}

	@Override
	public void init(double targetEpsilon) {

		status = Status.Estimating;
		
		this.targetEpsilon = targetEpsilon;

		values = new Values();

		energiedGMECEnumerated = false;
		scoredGMECPStarEnumerated = false;
		scoredGMECQPrimeEnumerated = false;

		// compute p*: boltzmann-weight the scores for all pruned conformations
		if(!computeGMECRatio) {
			ConfSearch ptree = confSearchFactory.make(emat, invmat);
			if(ptree instanceof ConfAStarTree) ((ConfAStarTree)ptree).stopProgress();
			values.pstar = calcWeightSumUpperBound(ptree);
		}

		// make the search tree for computing q*
		ConfSearch tree = confSearchFactory.make(emat, pmat);
		if(tree instanceof ConfAStarTree) ((ConfAStarTree)tree).stopProgress();
		ConfSearch.Splitter confsSplitter = new ConfSearch.Splitter(tree);
		scoreConfs = confsSplitter.makeStream();
		energyConfs = confsSplitter.makeStream();
		numConfsEvaluated = 0;
		numConfsToScore = tree.getNumConformations();		
		qprimeUnevaluated = BigDecimal.ZERO;
		qprimeUnscored = BigDecimal.ZERO;

		qstarScoreWeights = BigDecimal.ZERO;

		// treat the partition function state as if we have already processed 
		// the minGMEC
		if(minGMECConfs != null) {
			for(EnergiedConf conf : minGMECConfs) {
				values.qstar = values.qstar.add(boltzmann.calc(conf.getEnergy()));
				if (confListener != null) {
					confListener.onConf(conf);
				}
			}
			
			if(computeGMECRatio) {
				status = Status.Estimated;
				return;
			}

			numConfsEvaluated += minGMECConfs.size();
			numConfsToScore = numConfsToScore.subtract(BigInteger.valueOf(minGMECConfs.size()));

			//start up qprime before evaluating confs, because we would immediately
			//hit target epsilon otherwise
			values.qprime = updateQprime(null);
		}

		stopwatch = new Stopwatch().start();
	}

	protected BigDecimal calcWeightSumUpperBound(ConfSearch tree) {

		BigDecimal sum = BigDecimal.ZERO;
		BigDecimal boundOnAll = BigDecimal.ZERO;

		BigInteger numConfsRemaining = tree.getNumConformations();

		while (true) {

			// get the next conf
			ScoredConf conf = tree.nextConf();
			if (conf == null) {
				break;
			}

			/*
			//minGMECConfs are guaranteed not to be here, so no need to check here
			if(minGMECAssignments != null && !scoredGMECPStarEnumerated && contains(conf.getAssignments())) {
				scoredGMECPStarEnumerated = true;
				continue;
			}
			*/

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
			double tightness = boundOnRemaining.divide(boundOnAll, RoundingMode.HALF_UP).doubleValue();
			if (tightness <= 0.01) {
				break;
			}
		}

		return boundOnAll;
	}

	@Override
	protected BigDecimal updateQprime(EnergiedConf econf) {

		Stopwatch stopwatch = computeMaxNumConfs ? new Stopwatch().start() : null;

		// look through the conf tree to get conf scores
		// (which should be lower bounds on the conf energy)
		while (true) {

			// read a conf from the tree
			ScoredConf conf = scoreConfs.next();

			if (conf == null) {
				qprimeUnscored = BigDecimal.ZERO;
				break;
			}

			if(minGMECAssignments != null && !scoredGMECQPrimeEnumerated && contains(conf.getAssignments())) {
				scoredGMECQPrimeEnumerated = true;
				continue;
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
			double tightness = qprimeUnscored.divide(qprimeUnevaluated.add(qprimeUnscored), RoundingMode.HALF_UP).doubleValue();
			if (tightness <= 0.01) {
				break;
			}

			//reaching desired tightness can take a really long time
			//in a design, so ignore tightness if we are only interested
			//in computing the partition function to MaxNumConfs
			if(stopwatch != null && stopwatch.getTimeH() >= QPRIME_TIMEOUT_HRS) break;
		}

		if(econf != null) {
			qprimeUnevaluated = qprimeUnevaluated.subtract(boltzmann.calc(econf.getScore()));
		}
		return qprimeUnevaluated.add(qprimeUnscored);
	}
	
	boolean contains(int[] other) {
		ArrayList<Integer> assignment = new ArrayList<>();
		for(int rotamer : other) assignment.add(rotamer);
		return minGMECAssignments == null ? false : minGMECAssignments.contains(assignment);
	}

	protected ScoredConf getScoredConf() {
		ScoredConf conf;

		//don't double count the mingmec
		while (true) {

			conf = energyConfs.next();

			if (conf == null) {
				if(status != Status.Estimated) status = Status.NotEnoughConformations;
				return null;
			}

			//skip mingmec if it has already been enumerated
			if(minGMECAssignments != null && !energiedGMECEnumerated && contains(conf.getAssignments())) {
				energiedGMECEnumerated = true;
				continue;
			}

			if (boltzmann.calc(conf.getScore()).compareTo(BigDecimal.ZERO) == 0) {
				if(status != Status.Estimated) status = Status.NotEnoughFiniteEnergies;
				return null;
			}

			break;
		}

		return conf;
	}

	protected ArrayList<ScoredConf> getScoredConfs(int numRequested) {
		ScoredConf conf;
		if(scoredConfs == null) scoredConfs = new ArrayList<>(Collections.nCopies(numRequested, null));
		if(energiedConfs == null) energiedConfs = new ArrayList<>();
		int numScored = 0;

		while(numScored < numRequested) {
			if((conf = getScoredConf()) == null)
				break;
			scoredConfs.set(numScored++, conf);
		}

		//trim if numScored < numRequested
		if(numScored < numRequested) {
			scoredConfs.subList(numScored, scoredConfs.size()).clear();
		}

		return numScored > 0 ? scoredConfs : null;
	}

	protected void handleEnergiedConf(EnergiedConf econf) {
		if(status != Status.Estimated) {
			// get the boltzmann weight
			BigDecimal energyWeight = boltzmann.calc(econf.getEnergy());

			// update pfunc state
			numConfsEvaluated++;
			values.qstar = values.qstar.add(energyWeight);
			values.qprime = updateQprime(econf);

			// report progress if needed
			if (isReportingProgress && numConfsEvaluated % ecalc.getParallelism() == 0) {
				confOutput(econf);
			}

			// report confs if needed
			if (confListener != null) {
				confListener.onConf(econf);
			}

			// update status if needed
			double effectiveEpsilon = values.getEffectiveEpsilon();
			if(Double.isNaN(effectiveEpsilon)) {
				status = Status.NotEnoughFiniteEnergies;
			}
			else if (effectiveEpsilon <= targetEpsilon) {
				status = Status.Estimated;
				if (isReportingProgress) confOutput(econf);//just to let the user know we reached epsilon
			}
		}
	}

	@Override
	public void compute(long maxNumConfs) {

		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		long stopAtConf = numConfsEvaluated + maxNumConfs;
		while (true) {

			// get a conf from the tree
			// lock though to keep from racing the listener thread on the conf tree
			ScoredConf conf = null;
			ArrayList<ScoredConf> scoredConfs = null;

			synchronized (this) {

				// did we win?
				boolean hitEpsilonTarget = values.getEffectiveEpsilon() <= targetEpsilon;
				if (hitEpsilonTarget) {
					status = Status.Estimated;
					break;
				}

				// should we keep going?
				else if (!status.canContinue() || numConfsEvaluated >= stopAtConf) {
					break;
				}

				// get a list of scored confs
				else if (SYNCHRONIZED_MINIMIZATION) {
					scoredConfs = getScoredConfs(ecalc.getParallelism());
				}

				// or a single scored conf
				else {
					conf = getScoredConf();
				}
			}

			if (!status.canContinue()) {
				// we hit a failure condition and need to stop
				// but wait for current async minimizations to finish in case that pushes over the epsilon target
				// NOTE: don't wait for finish inside the sync, since that will definitely deadlock
				ecalc.waitForFinish();
				if(scoredConfs == null) {
					continue;
				}
			}

			// now handle scored confs
			// wait for all confs to finish minimizing, then process confs in order of score
			if (SYNCHRONIZED_MINIMIZATION) {
				
				for (ScoredConf sconf : scoredConfs) {
					ecalc.calcEnergyAsync(sconf, (EnergiedConf econf) -> {
						synchronized (this) {
							energiedConfs.add(econf);
						}
					});
				}

				while(energiedConfs.size() != scoredConfs.size()) {
					try { Thread.sleep(1); } catch (Exception e) {
						e.printStackTrace();
						System.exit(1);
					}
					ecalc.waitForFinish();
				}

				// sort energied confs by score
				Collections.sort(energiedConfs, new Comparator<EnergiedConf>() {
					@Override
					public int compare(EnergiedConf o1, EnergiedConf o2) {
						return o1.getScore() <= o2.getScore() ? -1 : 1;
					}
				});

				for (EnergiedConf econf : energiedConfs) {
					handleEnergiedConf(econf);
				}
				
				energiedConfs.clear();
			}

			// or do the energy calculation asynchronously
			else {			
				ecalc.calcEnergyAsync(conf, (EnergiedConf econf) -> {
					// energy calculation done
					// this is (potentially) running on a task executor listener thread
					// so lock to keep from racing the main thread
					synchronized (this) {
						handleEnergiedConf(econf);
					}
				});
			}
		}

		// wait for any remaining async minimizations to finish
		ecalc.waitForFinish();
	}

	void confOutput(ScoredConf conf) {
		MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
		double confVal = conf instanceof EnergiedConf ? ((EnergiedConf)conf).getEnergy() : conf.getScore();
		System.out.println(String.format("conf: %4d, energy: %.6f, q*: %12e, q': %12e, p*: %12e, range: [%12e, %12e], epsilon: %.6f, time: %10s, heapMem: %.0f%%",
				numConfsEvaluated, confVal, values.qstar, values.qprime, values.pstar, values.qstar, getUpperBound(), 
				values.getEffectiveEpsilon(),
				stopwatch.getTime(2),
				100f*heapMem.getUsed()/heapMem.getMax()
				));
	}

	public void setStatus(Status val) {
		status = val;
	}

	public void setValues(Values val) {
		this.values = val;
	}

	public void cleanup() {
		scoreConfs = null;
		energyConfs = null;

		scoredConfs = null;
		energiedConfs = null;

		confSearchFactory = null;
		confListener = null;
		
		emat = null;
		pmat = null;
		
		minGMECConfs = null;
		minGMECAssignments = null;
	}

	public void setNumConfsEvaluated(long val) {
		this.numConfsEvaluated = val;
	}

	public long getNumConfsEvaluated() {
		return this.numConfsEvaluated;
	}
	
	private long getNumPStarConfs() {
		ConfSearch ptree = confSearchFactory.make(emat, invmat);
		return ptree.getNumConformations().longValue();
	}
	
	private BigDecimal ExpE0() {
		ConfSearch ptree = confSearchFactory.make(emat, invmat);
		if(ptree instanceof ConfAStarTree) ((ConfAStarTree)ptree).stopProgress();
		
		ScoredConf conf = ptree.nextConf();
		
		if (conf == null) {
			return BigDecimal.ONE;
		}
		
		// compute the boltzmann weight for this conf
		BigDecimal weight = boltzmann.calc(conf.getScore());
		return weight;
	}
	
	private double rho() {
		return this.targetEpsilon/(1.0-this.targetEpsilon);
	}
	
	public long computeAdditionalConfs() {
		BigDecimal k = new BigDecimal(getNumPStarConfs());
		BigDecimal quotient = new BigDecimal(this.rho()).multiply(values.qstar).divide(this.ExpE0());
		BigDecimal ans = k.subtract(quotient);
		return ans.longValue()+1;
	}

	public void setComputeGMECRatio(boolean val) {
		this.computeGMECRatio = val;
	}

	public ConfListener getConfListener() {
		return this.confListener;
	}

	public void setMinGMECConfs(List<EnergiedConf> val) {
		this.minGMECConfs = val;
		if(val == null) return;
		
		this.minGMECAssignments = new HashSet<>();
		for(EnergiedConf conf : val) {
			ArrayList<Integer> assignment = new ArrayList<>();
			for(Integer rotamer : conf.getAssignments()) assignment.add(rotamer);
			assignment.trimToSize();
			this.minGMECAssignments.add(assignment);
		}
	}
	
	public List<EnergiedConf> getMinGMECConfs() {
		return this.minGMECConfs;
	}

	public BoltzmannCalculator getBoltzmannCalculator() {
		return this.boltzmann;
	}

	public void setComputeMaxNumConfs(boolean val) {
		this.computeMaxNumConfs = val;
	}
	
	public boolean getComputeMaxNumConfs() {
		return this.computeMaxNumConfs;
	}
	
	public BigDecimal getUpperBound() {
		return values.qstar.add(values.qprime).add(values.pstar);
	}
	
	public BigDecimal getLowerBound() {
		return values.qstar;
	}
	
	public String toString() {
		return String.format("q*: %12e, q': %12e, p*: %12e, range[%12e, %12e]", 
				values.qstar, values.qprime, values.pstar, values.qstar, 
				getUpperBound());
	}
}
