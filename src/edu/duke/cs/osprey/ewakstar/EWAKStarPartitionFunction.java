package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEPfunc;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.function.Function;

/** based on PartitionFunction.java, author: lowegard **/

public interface EWAKStarPartitionFunction {
	
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
				buf.append(" (log10)");
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

	void init(double targetEnergy, double targetEpsilon, BigInteger numConfsBeforePruning);
	/**
	 * Initializes the partition function for calculation.
	 * @param targetEpsilon The accuracy with which to estimate the partition function.
	 *
	 * @param stabilityThreshold If the partition function upper bound value falls below
	 *                           this threshold, the sequence is considered unstable.
	 * @param targetEnergy The energy window with which to estimate the partition function.
	 */
	default void init(double targetEnergy, double targetEpsilon, BigDecimal stabilityThreshold) {
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


	public static interface WithConfTable extends EWAKStarPartitionFunction {

		void setConfTable(ConfDB.ConfTable table);

		public static void setOrThrow(EWAKStarPartitionFunction pfunc, ConfDB.ConfTable table) {
			if (pfunc instanceof EWAKStarPartitionFunction.WithConfTable) {
				((EWAKStarPartitionFunction.WithConfTable)pfunc).setConfTable(table);
			} else {
				throw new EWAKStarPartitionFunction.WithConfTable.UnsupportedException(pfunc);
			}
		}

		public static class UnsupportedException extends RuntimeException {
			public UnsupportedException(EWAKStarPartitionFunction pfunc) {
				super("This partition function implementation (" + pfunc.getClass().getSimpleName() + ") doesn't support conformation database tables");
			}
		}
	}

	/**
	 * Factory method to make the best pfunc calculator based on the conf ecalc
	 */
	public static EWAKStarPartitionFunction makeBestFor(ConfSearch confSearch, ConfEnergyCalculator confEcalc) {

		// algorithms based on energy bounds can use the GD calculator, it's the most recent pfunc calculator
		return new EWAKStarGradientDescentPfunc(confSearch, confEcalc);
	}

	public static interface WithExternalMemory extends EWAKStarPartitionFunction {

		void setUseExternalMemory(boolean val, RCs rcs);

		public static void setOrThrow(EWAKStarPartitionFunction pfunc, boolean val, RCs rcs) {
			if (pfunc instanceof EWAKStarPartitionFunction.WithExternalMemory) {
				((EWAKStarPartitionFunction.WithExternalMemory)pfunc).setUseExternalMemory(val, rcs);
			} else {
				throw new EWAKStarPartitionFunction.WithExternalMemory.UnsupportedException(pfunc);
			}
		}

		public static class UnsupportedException extends RuntimeException {
			public UnsupportedException(EWAKStarPartitionFunction pfunc) {
				super("This partition function implementation (" + pfunc.getClass().getSimpleName() + ") doesn't support external memory");
			}
		}
	}
}
