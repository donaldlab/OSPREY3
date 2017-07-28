package edu.duke.cs.osprey.kstar.pfunc;

import java.math.BigDecimal;
import java.math.RoundingMode;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

public interface PartitionFunction {
	
	public static enum Status {
		
		Estimating(true),
		Estimated(false),
		OutOfConformations(false),
		OutOfLowEnergies(false);
		
		private boolean canContinue;
		
		private Status(boolean canContinue) {
			this.canContinue = canContinue;
		}
		
		public boolean canContinue() {
			return canContinue;
		}
	}
	
	public static class Values {
		
		public BigDecimal qstar; // pfunc value of all evaluated confs
		public BigDecimal qprime; // pfunc value of all unpruned, but unevaluated confs
		public BigDecimal pstar; // pfunc value of all pruned confs
		
		public Values() {
			qstar = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			qprime = BigDecimal.ZERO;
			pstar = BigDecimal.ZERO;
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

			BigDecimal s = qprime.add(pstar);
			BigDecimal qu = s.add(qstar);
			
			if (qu.compareTo(BigDecimal.ZERO) == 0) {
				// this is really bad... it should never happen
				// it means the upper bound on the pfunc value is zero
				// it probably means there are zero conformations in the tree?
				return Double.NaN;
			}
			
			return s.divide(qu, RoundingMode.HALF_UP).doubleValue();
		}

		public BigDecimal calcLowerBound() {
			return qstar;
		}

		public BigDecimal calcUpperBound() {
			return qstar.add(qprime).add(pstar);
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
			if (status == PartitionFunction.Status.Estimated) {
				return String.format("%e", values.qstar.doubleValue());
			}
			return status.toString();
		}
	}
	
	public static interface ConfListener {
		void onConf(ScoredConf conf);
	}
	
	void setReportProgress(boolean val);
	void setConfListener(ConfListener val);
	
	void init(double targetEpsilon);
	
	Status getStatus();
	Values getValues();
	int getParallelism();
	int getNumConfsEvaluated();
	
	void compute();
	void compute(int maxNumConfs);

	public default Result makeResult() {
		return new Result(getStatus(), getValues(), getNumConfsEvaluated());
	}
}
