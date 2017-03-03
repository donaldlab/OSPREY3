package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.ParallelConfPartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;

public class ParallelPartitionFunction extends ParallelConfPartitionFunction {

	protected PriorityQueue<EnergiedConf> econfs;
	protected int maxNumTopConfs;
	protected BigDecimal qstarScoreWeights;
	protected int numActiveThreads;

	public ParallelPartitionFunction(EnergyMatrix emat, PruningMatrix pmat, ConfSearchFactory confSearchFactory,
			Async ecalc) {
		super(emat, pmat, confSearchFactory, ecalc);
		qstarScoreWeights = null;
		econfs = null;
	}

	protected void writeTopConfs(int state, MSSearchProblem search) {
		if(econfs==null || econfs.size()==0) return;
		String seq = search.settings.getFormattedSequence();
		if(isReportingProgress) {
			System.out.println("Writing top "+ econfs.size()+" confs:");
		}
		seq = seq.replace(" ", ".");
		String dir = "topConfs"+File.separator+"State."+state+File.separator+seq;
		ObjectIO.makeDir(dir, false);
		for(int i=econfs.size()-1;i>-1;--i) {
			if(isReportingProgress)
				System.out.println(String.format("conf: %4d.pdf, energy: %.6f", i, econfs.peek().getEnergy()));
			String PDBFileName = dir+File.separator+i+".pdb";
			search.outputMinimizedStruct(econfs.poll().getAssignments(), PDBFileName);
		}
	}

	protected void saveEConf(EnergiedConf econf) {
		if(econfs.size() >= maxNumTopConfs) {
			if(econfs.peek().getEnergy() > econf.getEnergy()) econfs.poll();
			else return;
		}
		econfs.add(econf);
	}

	protected void saveEConfs(PriorityQueue<EnergiedConf> other) {
		if(econfs==null || other==null) return;
		while(other.size()>0) 
			saveEConf(other.poll());
	}

	@Override
	public void init(double targetEpsilon) {
		super.init(targetEpsilon);
		qstarScoreWeights = BigDecimal.ZERO;
		numActiveThreads = 0;
		maxNumTopConfs = 0;
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
						System.out.println(String.format("conf: %4d, energy: %.6f, q*: %12e, q': %12e, score diff: %12e, epsilon: %.6f, time: %10s, heapMem: %.0f%%",
								numConfsEvaluated, econf.getEnergy(), values.qstar, values.qprime, pdiff, values.getEffectiveEpsilon(),
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
