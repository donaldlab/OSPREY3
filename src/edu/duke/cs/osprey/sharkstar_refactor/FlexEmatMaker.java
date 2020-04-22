package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.io.File;

class FlexEmatMaker {

    public static EnergyMatrix makeEmat(ConfEnergyCalculator confECalc, String name, String cachePattern) {
        EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confECalc)
                //.setCacheFile(new File(name+".emat"))
                .setCacheFile(new File(cachePattern+"."+name+".emat"))
                .build()
                .calcEnergyMatrix();
        return emat;
    }

    public static ConfEnergyCalculator makeRigidConfEcalc(ConfEnergyCalculator sourceConfEcalc){
        EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(sourceConfEcalc.ecalc)
                .setIsMinimizing(false)
                .build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(sourceConfEcalc, ecalcRigid);
        return confEcalcRigid;
    }

    public static ConfEnergyCalculator makeMinimizeConfEcalc(SimpleConfSpace confSpace, Parallelism parallelism) {
        ForcefieldParams ffparams = new ForcefieldParams();

        // how should we compute energies of molecules?
        EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();
        // how should we define energies of conformations?
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                        .build()
                        .calcReferenceEnergies()
                )
                .build();
        return confEcalcMinimized;
    }
}
