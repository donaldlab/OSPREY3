package edu.duke.cs.osprey.kstar.pfunc;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.function.Function;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.kstar.KStarScore;
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
			BigDecimal x = MathTools.bigAdd(qstar, qprime, decimalPrecision);
			x = MathTools.bigAdd(x, pstar, decimalPrecision);
			return x;
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
	
	void init(double targetEpsilon);

	/**
	 * Initializes the partition function for calculation.
	 * @param targetEpsilon The accuracy with which to estimate the partition function.
	 * @param stabilityThreshold If the partition function upper bound value falls below
	 *                           this threshold, the sequence is considered unstable.
	 */
	default void init(double targetEpsilon, BigDecimal stabilityThreshold) {
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
	}
}
