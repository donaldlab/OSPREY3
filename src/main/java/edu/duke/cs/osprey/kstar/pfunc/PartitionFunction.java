/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.kstar.pfunc;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.time.Duration;
import java.util.function.Function;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;


public interface PartitionFunction {
	
	public static enum Status {
		
		Estimating(true),
		Estimated(false),
		OutOfConformations(false),
		OutOfLowEnergies(false),
		Unstable(false),
		Aborted(false);
		
		private boolean canContinue;
		
		private Status(boolean canContinue) {
			this.canContinue = canContinue;
		}
		
		public boolean canContinue() {
			return canContinue;
		}
	}

	public final MathContext decimalPrecision = new MathContext(64, RoundingMode.HALF_UP);
	
	public static class Values {
		
		public BigDecimal qstar; // pfunc value of all evaluated confs
		public BigDecimal qprime; // pfunc value of all unpruned, but unevaluated confs
		public BigDecimal pstar; // pfunc value of all pruned confs
		
		public Values() {
			qstar = BigDecimal.ZERO;
			qprime = BigDecimal.ZERO;
			pstar = BigDecimal.ZERO;

			// NOTE: q' should probably default to +inf
			// but our +inf implementaiton is a giant hack,
			// and will break existing code that's not expecting it
			// if your code wants the defualt +inf,
			// use the static factory method instead of this constructor
		}

		public static Values makeFullRange() {
			Values values = new Values();
			values.qprime = MathTools.BigPositiveInfinity;
			return values;
		}
		
		public double getEffectiveEpsilon() {

			// converting a single-sequence bound into an effective epsilon:
			
			// start with eqn 3 in K* paper:
			// q* = (1 - e)q
			
			// solve for e:
			// q* = q - qe
			// qe = q - q*
			// e = (q - q*)/q

			// by definition:
			// q := q* + p

			// upper bound on p based on pruned (p*) and un-enumerated (q') confs:
			// p <= s
			// s := q' + p*

			// upper bound on q (aka, the single-sequence bound):
			// q <= qu
			// qu := q* + s

			// substitute in upper bound on q to get "effective" epsilon, e*
			// NOTE: the BBK* paper calls this value delta
			// e* := (qu - q*)/qu

			// simplify:
			// e* = ((q* + s) - q*)/qu
			//    = s/qu

			BigDecimal s = MathTools.bigAdd(qprime, pstar, decimalPrecision);
			BigDecimal qu = MathTools.bigAdd(s, qstar, decimalPrecision);
			double delta = MathTools.bigDivide(s, qu, decimalPrecision).doubleValue();

			// sometimes the delta is ever so slightly below zero, so just take the abs
			delta = Math.abs(delta);

			return delta;
		}

		public BigDecimal calcLowerBound() {
			return qstar;
		}

		public BigDecimal calcUpperBound() {
			return new BigMath(decimalPrecision)
				.set(qstar)
				.add(qprime)
				.add(pstar)
				.get();
		}

		public MathTools.BigDecimalBounds calcBounds() {
			return new MathTools.BigDecimalBounds(
				calcLowerBound(),
				calcUpperBound()
			);
		}

		public double calcFreeEnergyLowerBound() {
			return new BoltzmannCalculator(PartitionFunction.decimalPrecision).freeEnergy(calcUpperBound());
		}

		public double calcFreeEnergyUpperBound() {
			return new BoltzmannCalculator(PartitionFunction.decimalPrecision).freeEnergy(calcLowerBound());
		}

		public MathTools.DoubleBounds calcFreeEnergyBounds(MathTools.DoubleBounds dest) {
			dest.lower = calcFreeEnergyLowerBound();
			dest.upper = calcFreeEnergyUpperBound();
			return dest;
		}

		public MathTools.DoubleBounds calcFreeEnergyBounds() {
			return new MathTools.DoubleBounds(
				calcFreeEnergyLowerBound(),
				calcFreeEnergyUpperBound()
			);
		}

		public double calcFreeEnergyLowerBoundPrecise() {
			return new BoltzmannCalculator(PartitionFunction.decimalPrecision).freeEnergyPrecise(calcUpperBound());
		}

		public double calcFreeEnergyUpperBoundPrecise() {
			return new BoltzmannCalculator(PartitionFunction.decimalPrecision).freeEnergyPrecise(calcLowerBound());
		}

		public MathTools.DoubleBounds calcFreeEnergyBoundsPrecise() {
			return new MathTools.DoubleBounds(
				calcFreeEnergyLowerBoundPrecise(),
				calcFreeEnergyUpperBoundPrecise()
			);
		}

		@Override
		public String toString() {
			return String.format("Z = %s", calcBounds());
		}
	}

	public static class Result {

		public final Status status;
		public final Values values;
		public final int numConfs;
		public final ConfListener confListener;

		public Result(Status status, Values values, int numConfs) {
			this(status, values, numConfs, null);
		}

		public Result(Status status, Values values, int numConfs, ConfListener listener) {
			this.status = status;
			this.values = values;
			this.numConfs = numConfs;
			this.confListener = listener;
		}

		@Override
		public String toString() {
			Function<String,String> trim = (s) -> {
				if (s.length() > 9) {
					return s.substring(0, 9);
				} else {
					return s;
				}
			};
			StringBuilder buf = new StringBuilder();
			buf.append(String.format("[%-9s,%9s]",
				trim.apply(KStarScore.scoreToLog10String(values.calcLowerBound())),
				trim.apply(KStarScore.scoreToLog10String(values.calcUpperBound()))
			));
			if (status == Status.Estimated) {
				buf.append(String.format(" %-26s", "(log10)"));
			} else {
				buf.append(String.format(" %-26s", "(log10," + status.name() + ")"));
			}
			return buf.toString();
		}

		public static Result makeAborted() {
			return new Result(Status.Aborted, Values.makeFullRange(), 0);
		}
	}
	
	public static interface ConfListener {
		void onConf(ScoredConf conf);
	}
	
	void setReportProgress(boolean val);
	void setConfListener(ConfListener val);

	/**
	 * Initializes the partition function for calculation.
	 * @param targetEpsilon The accuracy with which to estimate the partition function.
	 */
	void init(double targetEpsilon);

	/**
	 * Sets the stability threshold for this PartitionFunction, if supported
	 * @param stabilityThreshold If the partition function upper bound value falls below
	 *                           this threshold, the sequence is considered unstable.
	 */
	default void setStabilityThreshold(BigDecimal stabilityThreshold) {
		throw new UnsupportedOperationException(getClass().getName() + " does not yet support stability thresholds");
	}

	Status getStatus();
	Values getValues();
	int getParallelism();
	int getNumConfsEvaluated();

	void compute(int maxNumConfs);

	default void compute() {
		compute(Integer.MAX_VALUE);
	}

	/**
	 * Compute the partition function to the given epsilon,
	 * but stop earaly if the timeout was reached.
	 */
	default void compute(Duration timeout) {
		long startNs = System.nanoTime();
		long timeoutNs = timeout.toNanos();
		while (getStatus().canContinue) {
			// use the parallelism as the batch size
			// that's probably good, right? YOLO!
			int batchSize = getParallelism();
			compute(batchSize);
			long elapsedNs = System.nanoTime() - startNs;
			if (elapsedNs >= timeoutNs) {
				break;
			}
		}
	}

	public default Result makeResult() {
		return new Result(getStatus(), getValues(), getNumConfsEvaluated());
	}

	interface WithConfDB extends PartitionFunction {

		void setConfDB(ConfDB confDB, ConfDB.Key key);

		default void setConfDB(ConfDB confDB, String table) {
			setConfDB(confDB, new ConfDB.Key(table));
		}

		default void setConfDB(ConfDB confDB, Sequence seq) {
			setConfDB(confDB, new ConfDB.Key(seq));
		}

		/**
		 * Try to cast the pfunc to WithConfDB,
		 * but throw a nice error if the cast fails.
		 */
		static WithConfDB cast(PartitionFunction pfunc) {
			if (pfunc instanceof PartitionFunction.WithConfDB) {
				return (PartitionFunction.WithConfDB)pfunc;
			} else {
				throw new UnsupportedOperationException(
					"This partition function implementation (" + pfunc.getClass().getSimpleName() + ") doesn't support conformation databases"
				);
			}
		}
	}

	public static interface WithExternalMemory extends PartitionFunction {

		void setUseExternalMemory(boolean val, RCs rcs);

		public static void setOrThrow(PartitionFunction pfunc, boolean val, RCs rcs) {
			if (pfunc instanceof PartitionFunction.WithExternalMemory) {
				((PartitionFunction.WithExternalMemory)pfunc).setUseExternalMemory(val, rcs);
			} else {
				throw new PartitionFunction.WithExternalMemory.UnsupportedException(pfunc);
			}
		}

		public static class UnsupportedException extends RuntimeException {
			public UnsupportedException(PartitionFunction pfunc) {
				super("This partition function implementation (" + pfunc.getClass().getSimpleName() + ") doesn't support external memory");
			}
		}
	}

	/** Override to support task contexts, for contextual task executors */
	default void setInstanceId(int val) {
		// ignored by default
	}

	/** Override to support task contexts, for contextual task executors */
	default void putTaskContexts(TaskExecutor.ContextGroup contexts) {
		// ignored by default
	}
}
