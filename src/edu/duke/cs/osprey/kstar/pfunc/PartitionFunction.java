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
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.function.Function;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEPfunc;
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
			return MathTools.bigDivide(s, qu, decimalPrecision).doubleValue();
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
	}

	public static class Result {

		public final Status status;
		public final Values values;
		public final int numConfs;

		public Result(Status status, Values values, int numConfs) {
			this.status = status;
			this.values = values;
			this.numConfs = numConfs;
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
	 *
	 * Deprecated in favor of the two ConfSearch version, which allows computing pfuncs in constant memory
	 *
	 * @param confSearch The A* tree of conformations to enumerate (which may have been pruned)
	 * @param numConfsBeforePruning The total number of conformations in the conformation space for this search,
	 *                               including any conformations removed by pruned tuples.
	 * @param targetEpsilon The accuracy with which to estimate the partition function.
	 */
	@Deprecated
	void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon);

	/**
	 * Initializes the partition function for calculation using separate trees for computing upper bounds,
	 * and for lower bounds. Avoids the need to buffer conformations between the upper and lower bound calculations
	 * (either in internal or external memory), and allows running partition function calculations in constant memory.
	 * @param upperBoundConfs The A* tree of conformations to enumerate for computing upper bounds
	 * @param lowerBoundConfs The A* tree of conformations to enumerate for computing lower bounds
	 * @param numConfsBeforePruning The total number of conformations in the conformation space for this search,
	 *                               including any conformations removed by pruned tuples.
	 * @param targetEpsilon The accuracy with which to estimate the partition function.
	 */
	default void init(ConfSearch upperBoundConfs, ConfSearch lowerBoundConfs, BigInteger numConfsBeforePruning, double targetEpsilon) {
		throw new UnsupportedOperationException("This pfunc calculator (" + getClass().getSimpleName() + ") doesn't support two-tree initialization.");
	}

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

	public default Result makeResult() {
		return new Result(getStatus(), getValues(), getNumConfsEvaluated());
	}


	public static interface WithConfTable extends PartitionFunction {

		void setConfTable(ConfDB.ConfTable table);

		public static void setOrThrow(PartitionFunction pfunc, ConfDB.ConfTable table) {
			if (pfunc instanceof PartitionFunction.WithConfTable) {
				((PartitionFunction.WithConfTable)pfunc).setConfTable(table);
			} else {
				throw new PartitionFunction.WithConfTable.UnsupportedException(pfunc);
			}
		}

		public static class UnsupportedException extends RuntimeException {
			public UnsupportedException(PartitionFunction pfunc) {
				super("This partition function implementation (" + pfunc.getClass().getSimpleName() + ") doesn't support conformation database tables");
			}
		}
	}

	/**
	 * Factory method to make the best pfunc calculator based on the conf ecalc
	 */
	public static PartitionFunction makeBestFor(ConfEnergyCalculator confEcalc) {
		if (confEcalc instanceof LUTEConfEnergyCalculator) {
			// LUTE needs its own calculator, since it doesn't use energy bounds
			return new LUTEPfunc((LUTEConfEnergyCalculator)confEcalc);
		} else {
			// algorithms based on energy bounds can use the GD calculator, it's the most recent pfunc calculator
			return new GradientDescentPfunc(confEcalc);
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
}
