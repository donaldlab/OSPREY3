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

	public double sum() {
    	double sum = 0.0;
		for (int i=0; i<oneBody.length; i++) {
			sum += oneBody[i];
		}
		for (int i=0; i<pairwise.length; i++) {
			sum += pairwise[i];
		}
		return sum;
	}

	@Override
	public String toString() {
    	return toString(6, (energy) -> {
			if (energy == Double.POSITIVE_INFINITY) {
				return String.format("%6s", "inf");
			} else if (energy == Double.NEGATIVE_INFINITY) {
				return String.format("%6s", "-inf");
			} else if (energy == 1000) {
				return String.format("%6s", "1k");
			} else if (energy > 1000) {
				return String.format("%6s", ">1k");
			} else if (energy == -1000) {
				return String.format("%6s", "-1k");
			} else if (energy < -1000) {
				return String.format("%6s", "<-1k");
			} else {
				return String.format("%6.2f", energy);
			}
		});
	}

	public String toString(int cellWidth, int precision) {
    	return toString(cellWidth, (energy) -> {
			if (energy == Double.POSITIVE_INFINITY) {
				return String.format("%" + cellWidth + "s", "inf");
			} else if (energy == Double.NEGATIVE_INFINITY) {
				return String.format("%" + cellWidth + "s", "-inf");
			} else {
				return String.format("%" + cellWidth + "." + precision + "f", energy);
			}
		});
	}

	public String toStringScientific() {
		return toString(8, (energy) -> {
			if (energy == Double.POSITIVE_INFINITY) {
				return String.format("%6s", "inf");
			} else if (energy == Double.NEGATIVE_INFINITY) {
				return String.format("%6s", "-inf");
			} else {
				return String.format("%6.1e", energy);
			}
		});
	}
}
