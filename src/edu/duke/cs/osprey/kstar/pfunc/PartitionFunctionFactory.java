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
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.lute.LUTEPfunc;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound;
import edu.duke.cs.osprey.sharkstar.SHARKStarBound;

import java.math.BigInteger;
import java.util.HashMap;
import java.util.Map;

public class PartitionFunctionFactory {

    private Map<String, MultiSequenceSHARKStarBound> stateBounds = new HashMap<>();

    enum PartitionFunctionImpl {
        MARKStar,
        GradientDescent,
        LUTE,
        SHARKStar,
        MSSHARKStar,
    }

    private ConfEnergyCalculator confUpperBoundECalc;
    private ConfEnergyCalculator confEcalc;
    private Map<ConfEnergyCalculator, EnergyMatrix> emats = new HashMap<>();
    private SimpleConfSpace confSpace;
    private EnergyMatrix upperBoundEmat;
    private PartitionFunctionImpl pfuncImpl = PartitionFunctionImpl.GradientDescent;
    private UpdatingEnergyMatrix MARKStarEmat = null;
    private String state = "(undefined)";
    private SHARKStarBound preComputedFlex = null;
    private MultiSequenceSHARKStarBound preComputedMSFlex = null;
    private MultiSequenceSHARKStarBound fullMSBound = null;

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

    public void setUseSHARKStar(ConfEnergyCalculator rigidConfECalc, SHARKStarBound preComputedFlex){
        this.preComputedFlex = preComputedFlex;
        setUseSHARKStar(rigidConfECalc);
    }

    public void setUseSHARKStar(ConfEnergyCalculator rigidConfECalc) {
        this.confUpperBoundECalc = rigidConfECalc;
        this.pfuncImpl = PartitionFunctionImpl.SHARKStar;
    }

    public void setPrecomputedCorrections(UpdatingEnergyMatrix corrections){
        this.MARKStarEmat = corrections;
    }

    public void setUseMSSHARKStar(ConfEnergyCalculator rigidConfECalc, MultiSequenceSHARKStarBound preComputedFlex){
        this.preComputedMSFlex = preComputedFlex;
        setUseMSSHARKStar(rigidConfECalc);
    }

    public void setUseMSSHARKStar(ConfEnergyCalculator rigidConfECalc) {
        this.confUpperBoundECalc = rigidConfECalc;
        this.pfuncImpl = PartitionFunctionImpl.MSSHARKStar;
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
        return makePartitionFunctionFor(rcs, confSpaceSize, epsilon, null, null);
    }

    public PartitionFunction makePartitionFunctionFor(RCs rcs, BigInteger confSpaceSize, double epsilon, PruningMatrix pmat) {
        return makePartitionFunctionFor(rcs, confSpaceSize, epsilon, pmat,null);
    }

    public PartitionFunction makePartitionFunctionFor(RCs rcs, BigInteger confSpaceSize, double epsilon, Sequence seq) {
        return makePartitionFunctionFor(rcs, confSpaceSize, epsilon, null,seq);
    }

    public PartitionFunction makePartitionFunctionFor(double epsilon, Sequence seq) {
        if(pfuncImpl != PartitionFunctionImpl.MSSHARKStar)
            throw new UnsupportedOperationException("Can only generate sequence-specific partition funcitions with SHARK*");
        RCs seqRCs = seq.makeRCs(confSpace);
        return makePartitionFunctionFor(seqRCs, seqRCs.getNumConformations(), epsilon, seq);
    }

    public PartitionFunction makePartitionFunctionFor(RCs rcs, BigInteger confSpaceSize, double epsilon, PruningMatrix pmat, Sequence seq) {
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
            case SHARKStar:
                minimizingEmat = makeEmat(confEcalc, "minimizing");
                //if(MARKStarEmat == null)
                    //MARKStarEmat = new UpdatingEnergyMatrix(confSpace, minimizingEmat, confEcalc);
                SHARKStarBound SHARKStarBound = null;
                if(preComputedFlex == null) {
                    SHARKStarBound = new SHARKStarBound(confSpace, makeEmat(confUpperBoundECalc, "rigid"),
                            minimizingEmat, confEcalc, rcs, confEcalc.ecalc.parallelism);
                }
                else {
                    SHARKStarBound = new SHARKStarBound(confSpace, makeEmat(confUpperBoundECalc, "rigid"),
                            minimizingEmat, confEcalc, rcs, confEcalc.ecalc.parallelism, preComputedFlex);
                }
                if (MARKStarEmat != null)
                    SHARKStarBound.mergeCorrections(MARKStarEmat);
                SHARKStarBound.init(epsilon);
                pfunc = SHARKStarBound;
                break;
            case MSSHARKStar:
                String RCString = rcs.toString();
                if(!stateBounds.containsKey(RCString)) {
                    MultiSequenceSHARKStarBound stateMSBound;
                    minimizingEmat = makeEmat(confEcalc, "minimizing");
                    if(MARKStarEmat == null)
                        MARKStarEmat = new UpdatingEnergyMatrix(confSpace, minimizingEmat, confEcalc);
                    if (preComputedMSFlex == null) {
                        stateMSBound = new MultiSequenceSHARKStarBound(confSpace, makeEmat(confUpperBoundECalc, "rigid"),
                                minimizingEmat, confEcalc, rcs, confEcalc.ecalc.parallelism);
                    } else {
                        stateMSBound = new MultiSequenceSHARKStarBound(confSpace, makeEmat(confUpperBoundECalc, "rigid"),
                                minimizingEmat, confEcalc, rcs, confEcalc.ecalc.parallelism, preComputedMSFlex);
                    }
                    stateMSBound.setCorrections(MARKStarEmat);
                    stateMSBound.init(epsilon);
                    stateBounds.put(RCString, stateMSBound);
                }
                pfunc = stateBounds.get(RCString);
                if(seq != null)
                    pfunc = stateBounds.get(RCString).getPartitionFunctionForSequence(seq);
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
