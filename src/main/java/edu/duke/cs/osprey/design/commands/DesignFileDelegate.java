package edu.duke.cs.osprey.design.commands;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.converters.FileConverter;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.design.models.MoleculeDto;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.PDBIO;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaDeviceProp;

import java.io.File;
import java.util.stream.Collectors;

import static jcuda.runtime.JCuda.cudaGetDeviceCount;
import static jcuda.runtime.JCuda.cudaGetDeviceProperties;

public class DesignFileDelegate {

    @Parameter(description = "Path to design file.", names = {"--design", "-d"}, converter = FileConverter.class, validateWith=FileExistsValidation.class)
    public File design;

    @Parameter(description = "Prints this help information", names={"--help", "-h"}, help = true)
    boolean help;

    @Parameter(description = "The approximation accuracy. Z* = (1 - epsilon)Z", names={"--epsilon", "-e"})
    double epsilon;

    @Parameter(names = "--verify-design", description = "Verifies input parameters, but does not run the design")
    boolean verifyInput;

    @Parameter(names = {"--cuda", "-c"})
    boolean useCuda;

    int getNumGpu() {
        var deviceCount = new int[]{0};
        cudaGetDeviceCount(deviceCount);
        return deviceCount[0];
    }

    int getNumberGpuMultiprocessors(int deviceIdx) {
        var deviceProps = new cudaDeviceProp();
        cudaGetDeviceProperties(deviceProps, deviceIdx);
        return deviceProps.multiProcessorCount;
    }

    Parallelism getParallelism() {
        if (useCuda) {
            JCuda.setExceptionsEnabled(true);
            var numGpus = getNumGpu();
            var numMultiprocessors = getNumberGpuMultiprocessors(0);
            return new Parallelism(Runtime.getRuntime().availableProcessors(), numGpus, numMultiprocessors);
        }

        return new Parallelism(Runtime.getRuntime().availableProcessors(), 0, 0);
    }

    SimpleConfSpace createConfSpace(MoleculeDto spec, ForcefieldParams ffParams) {
        /* Reads a PDB file into a Molecule. */
        var molecule = PDBIO.read(spec.coordinates);

        /* Reads the templates, rotamers, and entropy info for a given forcefield */
        /*
            "/config/parm96a.dat",
            "/config/all_amino94.in",
            "/config/all_aminont94.in",
            "/config/all_aminoct94.in",
            "/config/all_nuc94_and_gr.in",
         */
        var extraTemplates = spec.extraTemplates;
        var extraTemplateCoords = spec.extraTemplatesCoordinates;
        var additionalRotamers = spec.additionalRotamers;

        var templateLibrary = new ResidueTemplateLibrary.Builder(ffParams.forcefld)
                .addTemplates(extraTemplates)
                .addTemplateCoords(extraTemplateCoords)
                .addRotamers(additionalRotamers)
                .build();

        /* Strands combine a Molecule with design flexibility and templates */
        var strandBuilder = new Strand.Builder(molecule)
                .setTemplateLibrary(templateLibrary);

        /* Add in flexibility and mutability parameters */
        for (var mod : spec.residueModifiers) {
            var identifier = mod.identity.positionIdentifier();
            strandBuilder.setResidueMutability(
                    identifier,
                    mod.mutability,
                    mod.flexibility.includeStructureRotamer,
                    mod.flexibility.useContinuous
            );
        }
        var protein = strandBuilder.build();

        // Validate that design intentions match input structure
        var identities = spec.residueModifiers.stream().map(m -> m.identity).collect(Collectors.toUnmodifiableList());
        for(var identity : identities) {
            var flex = protein.flexibility.get(identity.positionIdentifier());
            if (!flex.wildType.substring(0, 2).equalsIgnoreCase(identity.aminoAcidType.substring(0, 2))) {
                throw new RuntimeException(String.format("Design parameter thought residue %s was %s, but in structure is %s", identity.positionIdentifier(), identity.aminoAcidType, flex.wildType));
            }
        }

        /* Maintains flexibility information with the molecule, and can use that to make new molecules */

        return new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();
    }

}
