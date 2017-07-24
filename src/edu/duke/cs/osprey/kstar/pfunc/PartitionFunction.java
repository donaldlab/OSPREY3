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
			
			// start with eqn 3 in K* paper
			// q* = (1 - e)q
			
			// solve for e:
			// q* = q - qe
			// qe = q - q*
			// e = (q - q*)/q
			
			// q := q* + q' + p*
			
			// so
			// e = (q* + q' + p* - q*)/(q* + q' + p*)
			// e = (q' + p*)/(q* + q' + p*)
			
			// for simplicity,
			// s := q' + p*
			
			// so
			// e = s/q
			
			BigDecimal s = qprime.add(pstar);
			BigDecimal q = s.add(qstar);
			
			if (q.compareTo(BigDecimal.ZERO) == 0) {
				// this is really bad... it should never happen
				// it probably means there are zero conformations in the tree
				return Double.NaN;
			}
			
			return s.divide(q, RoundingMode.HALF_UP).doubleValue();
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

		public static BigDecimal calcKStarScore(Result protein, Result ligand, Result complex) {
			if (protein.status == PartitionFunction.Status.Estimated
				&& ligand.status == PartitionFunction.Status.Estimated
				&& complex.status == PartitionFunction.Status.Estimated) {

				return complex.values.qstar
					.divide(protein.values.qstar, RoundingMode.HALF_UP)
					.divide(ligand.values.qstar, RoundingMode.HALF_UP);
			}
			return null;
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
