package edu.duke.cs.osprey.sparse;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.structure.Molecule;

import java.util.List;

public class SparseGraphGenerator {

    public static ResidueInteractions generateGraph(SimpleConfSpace confSpace, EnergyMatrix emat,
                                                    double energyCutoff, double distanceCutoff) {
        ResidueInteractions graph = new ResidueInteractions();
        for(SimpleConfSpace.Position pos1 : confSpace.positions) {
            for (SimpleConfSpace.Position pos2 : confSpace.positions) {
                if (pos1.equals(pos2))
                    continue;
                double minDistance = Double.MAX_VALUE;
                double minEnergy = Double.MAX_VALUE;
                double maxEnergy = Double.MIN_VALUE;
                for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
                    for (SimpleConfSpace.ResidueConf rc2 : pos2.resConfs) {
                        RCTuple pair = new RCTuple(pos1.index, rc1.index, pos2.index, rc2.index);
                        double energy = emat.getInternalEnergy(pair);
                        Molecule assignedMolecule = confSpace.makeDiscreteMolecule(pair);
                        double distance = assignedMolecule.getResByPDBResNumber(pos1.resNum)
                                .distanceTo(assignedMolecule.getResByPDBResNumber(pos2.resNum));
                        maxEnergy = Math.max(energy, maxEnergy);
                        minEnergy = Math.min(energy, minEnergy);
                        minDistance = Math.min(distance, minDistance);
                    }
                }
                if(maxEnergy - minEnergy > energyCutoff && minDistance < distanceCutoff)
                    graph.add(new ResidueInteractions.Pair(pos1.resNum, pos2.resNum,
                            ResidueInteractions.Pair.IdentityWeight, ResidueInteractions.Pair.IdentityOffset));
            }
            graph.addSingle(pos1.resNum);
        }
        return graph;
    }
}
