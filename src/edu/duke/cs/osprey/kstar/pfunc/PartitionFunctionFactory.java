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

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import java.io.File;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEPfunc;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.math.BigInteger;
import java.util.HashMap;
import java.util.Map;

public class PartitionFunctionFactory {

    enum PartitionFunctionImpl {
        MARKStar,
        GradientDescent,
        LUTE
    }

    private ConfEnergyCalculator confUpperBoundECalc;
    private ConfEnergyCalculator confEcalc;
    private Map<ConfEnergyCalculator, EnergyMatrix> emats = new HashMap<>();
    private SimpleConfSpace confSpace;
    private EnergyMatrix upperBoundEmat;
    private PartitionFunctionImpl pfuncImpl = PartitionFunctionImpl.GradientDescent;
    private UpdatingEnergyMatrix MARKStarEmat = null;
    private String state = "(undefined)";

    public PartitionFunctionFactory(SimpleConfSpace confSpace, ConfEnergyCalculator confECalc, String state) {
        this.state = state;
        this.confSpace = confSpace;
        this.confEcalc = confECalc;
    }

    public void setUseMARKStar(ConfEnergyCalculator rigidConfECalc) {
        this.confUpperBoundECalc = rigidConfECalc;
        this.pfuncImpl = PartitionFunctionImpl.MARKStar;
    }

    public void setUseLUTE(ConfEnergyCalculator confECalc) {
        this.confEcalc = confECalc;
        this.pfuncImpl = PartitionFunctionImpl.LUTE;
    }

    public void setUseGradientDescent() {
        this.pfuncImpl = PartitionFunctionImpl.GradientDescent;
    }

    public ConfSearch makeConfSearch(EnergyMatrix emat, RCs rcs, PruningMatrix pmat) {
        if(pmat != null)
            rcs = new RCs(rcs, pmat);
        return new ConfAStarTree.Builder(emat, rcs)
                .setTraditional()
                .build();
    }

    public ConfSearch makeConfSearch(EnergyMatrix emat, RCs rcs) {
        return makeConfSearch(emat, rcs, null);
    }


    public PartitionFunction makePartitionFunctionFor(RCs rcs, BigInteger confSpaceSize, double epsilon) {
        return makePartitionFunctionFor(rcs, confSpaceSize, epsilon, null);
    }

    public PartitionFunction makePartitionFunctionFor(RCs rcs, BigInteger confSpaceSize, double epsilon, PruningMatrix pmat) {
        PartitionFunction pfunc = null;
        switch (pfuncImpl) {
            case GradientDescent:
                pfunc = new GradientDescentPfunc(confEcalc);
                ConfSearch AStarSearch = makeConfSearch(makeEmat(confEcalc), rcs);
                pfunc.init(AStarSearch, rcs.getNumConformations(), epsilon);
                break;
            case MARKStar:
                EnergyMatrix minimizingEmat = makeEmat(confEcalc, "minimizing");
                if(MARKStarEmat == null)
                    MARKStarEmat = new UpdatingEnergyMatrix(confSpace, minimizingEmat, confEcalc);
                MARKStarBound MARKStarBound = new MARKStarBoundFastQueues(confSpace, makeEmat(confUpperBoundECalc, "rigid"),
                        minimizingEmat, confEcalc, rcs, confEcalc.ecalc.parallelism);
                MARKStarBound.setCorrections(MARKStarEmat);
                MARKStarBound.init(epsilon);
                pfunc = MARKStarBound;
                break;
            case LUTE:
                pfunc = new LUTEPfunc((LUTEConfEnergyCalculator) confEcalc);
                pfunc.init(
                        makeConfSearch(makeEmat(confEcalc), rcs),
                        makeConfSearch(makeEmat(confEcalc), rcs),
                        rcs.getNumConformations(),
                        epsilon
                );
                break;
        }
        return pfunc;

    }

    private EnergyMatrix makeEmat(ConfEnergyCalculator confECalc) {
        return makeEmat(confECalc, "default");
    }

    private EnergyMatrix makeEmat(ConfEnergyCalculator confECalc, String name) {
        if(!emats.containsKey(confECalc)) {
            System.out.println("Making energy matrix for "+confECalc);
            EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confECalc)
                    .setCacheFile(new File(state+"."+name+".emat"))
                    .build()
                    .calcEnergyMatrix();
            emats.put(confECalc, emat);
        }
        return emats.get(confECalc);
    }
}
