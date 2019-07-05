package edu.duke.cs.osprey.design;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.design.models.AminoAcid;
import edu.duke.cs.osprey.design.models.StabilityDesign;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.io.IOException;
import java.util.stream.Collectors;

public class PartitionFunctionDesign {

    public static void main(String[] args) throws IOException {
        var designFile = new File(args[0]);
        var stabilityDesign = StabilityDesign.parse(designFile);
        var molecule = PDBIO.read(stabilityDesign.molecule);

        var ffParams = new ForcefieldParams();
        var templateLibrary = new ResidueTemplateLibrary.Builder(ffParams.forcefld)
                .build();

        var protein = new Strand.Builder(molecule)
                .setTemplateLibrary(templateLibrary)
                .build();

        stabilityDesign.residueModifiers.forEach(mod -> {
            var residue = protein.flexibility.get(mod.identity.positionIdentifier());
            residue.addWildTypeRotamers = mod.flexibility.includeStructureRotamer;
            var toMutations = mod.mutable.stream().map(AminoAcid::toValue).collect(Collectors.toUnmodifiableList());
            if (!toMutations.isEmpty()) {
                residue.setLibraryRotamers(toMutations);
            }
            residue.setContinuous();
        });

        var confSpace = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .setShellDistance(4)
                .build();

        var parallelism = new Parallelism(Runtime.getRuntime().availableProcessors(), 0, 0);

        var energyCalculator = new EnergyCalculator.Builder(confSpace, ffParams)
                .setParallelism(parallelism)
                .build();

        var confEnergyCalculator = new ConfEnergyCalculator.Builder(confSpace, energyCalculator)
                .build();

        var partitionFnBuilder = new PartitionFunctionFactory(confSpace, confEnergyCalculator, "default");
        partitionFnBuilder.setUseGradientDescent();

        var rcs = new RCs(confSpace);
        var partFn = partitionFnBuilder.makePartitionFunctionFor(rcs, null, stabilityDesign.epsilon);
        partFn.compute();
        var evaluated = partFn.getNumConfsEvaluated();
        var i = 0;
    }
}
