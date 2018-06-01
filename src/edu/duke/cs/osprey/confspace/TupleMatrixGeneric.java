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

public class TupleMatrixGeneric<T> extends AbstractTupleMatrix<T> {

    private static final long serialVersionUID = 845854459137269739L;

    private ArrayList<T> oneBody; // indices: res1, RC1
    private ArrayList<T> pairwise; // indices: res1, res2, RC1, RC2 where res1>res2

    protected TupleMatrixGeneric() {
        // do nothing
        // so subclasses can avoid allocation if they want
        // TODO: make a boolean specialization, and get rid of this constructor
        // only the pruning matrix uses this functionality
    }

    public TupleMatrixGeneric(SimpleConfSpace confSpace) {
        this(confSpace.getNumPos(), confSpace.getNumResConfsByPos(), 0.0, null);
    }

    public TupleMatrixGeneric(ConfSpace cSpace, double pruningInterval, T defaultHigherInteraction) {
        super(cSpace, pruningInterval, defaultHigherInteraction);
    }

    public TupleMatrixGeneric(int numPos, int[] numConfAtPos, double pruningInterval, T defaultHigherInteraction){
        super(numPos, numConfAtPos, pruningInterval, defaultHigherInteraction);
    }

    @Override
    protected void allocate(int numOneBody, int numPairwise) {
        oneBody = new ArrayList<>(numOneBody);
        for (int i=0; i<numOneBody; i++) {
            oneBody.add(null);
        }
        pairwise = new ArrayList<>(numPairwise);
        for (int i=0; i<numPairwise; i++) {
            pairwise.add(null);
        }
    }

    @Override
    public T getOneBody(int res, int conf) {
        return oneBody.get(getOneBodyIndex(res, conf));
    }

    @Override
    public void setOneBody(int res, int conf, T val) {
        oneBody.set(getOneBodyIndex(res, conf), val);
    }

    @Override
    public void setOneBody(int res, ArrayList<T> val) {
        int n = getNumConfAtPos(res);
        for (int i=0; i<n; i++) {
            oneBody.set(getOneBodyIndex(res, i), val.get(i));
        }
    }

    @Override
    public T getPairwise(int res1, int conf1, int res2, int conf2) {
        return pairwise.get(getPairwiseIndex(res1, conf1, res2, conf2));
    }

    @Override
    public void setPairwise(int res1, int conf1, int res2, int conf2, T val) {
        pairwise.set(getPairwiseIndex(res1, conf1, res2, conf2), val);
    }

    @Override
    public void setPairwise(int res1, int res2, ArrayList<ArrayList<T>> val) {
        int n1 = getNumConfAtPos(res1);
        int n2 = getNumConfAtPos(res2);
        for (int i1=0; i1<n1; i1++) {
            for (int i2=0; i2<n2; i2++) {
                pairwise.set(getPairwiseIndex(res1, i1, res2, i2), val.get(i1).get(i2));
            }
        }
    }
}