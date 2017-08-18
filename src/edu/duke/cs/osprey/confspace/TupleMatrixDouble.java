/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

import org.apache.commons.collections4.iterators.ArrayIterator;

public class TupleMatrixDouble extends AbstractTupleMatrix<Double> {
	
	private static final long serialVersionUID = -1286639255089978027L;
	
    //note: tuples are sets not ordered pairs, i.e. E(i_r,j_s) = E(j_s,i_r), and pruning (i_r,j_s) means pruning (j_s,i_r)
	private double[] oneBody; // indices: res1, RC1
	private double[] pairwise; // indices: res1, res2, RC1, RC2 where res1>res2
    
    public TupleMatrixDouble(ConfSpace cSpace, double pruningInterval, double defaultHigherInteraction) {
    	super(cSpace, pruningInterval, defaultHigherInteraction);
    }
    
    public TupleMatrixDouble(SimpleConfSpace confSpace, double pruningInterval, double defaultHigherInteraction) {
        super(confSpace, pruningInterval, defaultHigherInteraction);
    }
    
    public TupleMatrixDouble(int numPos, int[] numAllowedAtPos, double pruningInterval, double defaultHigherInteraction) {
    	super(numPos, numAllowedAtPos, pruningInterval, defaultHigherInteraction);
    }
    
    public TupleMatrixDouble(TupleMatrixDouble other) {
    	super(other);
    	this.oneBody = other.oneBody.clone();
    	this.pairwise = other.pairwise.clone();
    }
    
    @Override
    protected void allocate(int numOneBody, int numPairwise) {
        oneBody = new double[numOneBody];
        pairwise = new double[numPairwise];
    }
    
    @Override
    public Double getOneBody(int res, int conf) {
    	return oneBody[getOneBodyIndex(res, conf)];
    }
    
    @Override
    public void setOneBody(int res, int conf, Double val) {
    	oneBody[getOneBodyIndex(res, conf)] = val;
    }
    
    @Override
    public void setOneBody(int res, ArrayList<Double> val) {
    	int n = getNumConfAtPos(res);
    	for (int i=0; i<n; i++) {
    		oneBody[getOneBodyIndex(res, i)] = val.get(i);
    	}
    }
    
    @Override
    public Double getPairwise(int res1, int conf1, int res2, int conf2) {
    	return pairwise[getPairwiseIndex(res1, conf1, res2, conf2)];
    }
    
    @Override
    public void setPairwise(int res1, int conf1, int res2, int conf2, Double val) {
    	pairwise[getPairwiseIndex(res1, conf1, res2, conf2)] = val;
    }
    
    @Override
    public void setPairwise(int res1, int res2, ArrayList<ArrayList<Double>> val) {
    	int n1 = getNumConfAtPos(res1);
    	int n2 = getNumConfAtPos(res2);
    	for (int i1=0; i1<n1; i1++) {
    		for (int i2=0; i2<n2; i2++) {
    			pairwise[getPairwiseIndex(res1, i1, res2, i2)] = val.get(i1).get(i2);
    		}
    	}
    }
    
    public void fill(double[] vals) {
    	ArrayIterator<Double> iter = new ArrayIterator<>(vals);
    	fill(iter);
    	assert (iter.hasNext() == false);
    }

	public void negate() {
		for (int i=0; i<oneBody.length; i++) {
			oneBody[i] = -oneBody[i];
		}
		for (int i=0; i<pairwise.length; i++) {
			pairwise[i] = -pairwise[i];
		}
	}

	@Override
	public String toString() {

    	StringBuilder buf = new StringBuilder();

    	final String spacer = "  ";
		BiConsumer<Integer,Integer> labelPrinter = (pos, rc) -> {
			buf.append(spacer);
			buf.append(String.format("%2d:%03d", pos, rc));
		};
		Consumer<Double> energyPrinter = (energy) -> {
			buf.append(spacer);
			if (energy == Double.POSITIVE_INFINITY) {
				buf.append(String.format("%6s", "inf"));
			} else if (energy == Double.NEGATIVE_INFINITY) {
				buf.append(String.format("%6s", "-inf"));
			} else if (energy >= 1000) {
				buf.append(String.format("%6s", ">1k"));
			} else {
				buf.append(String.format("%6.2f", energy));
			}
		};
		Runnable blankPrinter = () -> {
			buf.append(spacer);
			buf.append("      ");
		};

    	// singles
		buf.append("singles:\n");
		blankPrinter.run();
		for (int pos1=0; pos1<getNumPos(); pos1++) {
			int n1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<n1; rc1++) {
				labelPrinter.accept(pos1, rc1);
			}
		}
		buf.append("\n");
		blankPrinter.run();
		for (int pos1=0; pos1<getNumPos(); pos1++) {
			int n1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<n1; rc1++) {
				energyPrinter.accept(getOneBody(pos1, rc1));
			}
		}
		buf.append("\n");

		// pairs
		buf.append("pairs:\n");
		blankPrinter.run();
		for (int pos1=0; pos1<getNumPos(); pos1++) {
			int n1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<n1; rc1++) {
				labelPrinter.accept(pos1, rc1);
			}
		}
		buf.append("\n");
		for (int pos1=1; pos1<getNumPos(); pos1++) {
			int n1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<n1; rc1++) {

				labelPrinter.accept(pos1, rc1);

				for (int pos2=0; pos2<pos1; pos2++) {
					int n2 = getNumConfAtPos(pos2);
					for (int rc2=0; rc2<n2; rc2++) {
						energyPrinter.accept(getPairwise(pos1, rc1, pos2, rc2));
					}
				}

				buf.append("\n");
			}
		}

		return buf.toString();
	}
}
