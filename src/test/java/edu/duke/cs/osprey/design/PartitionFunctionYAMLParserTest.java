package edu.duke.cs.osprey.design;

import edu.duke.cs.osprey.design.models.StabilityDesign;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertNotNull;

class PartitionFunctionYAMLParserTest {

    private StabilityDesign stabilityDesign;

    @BeforeEach
    void getTestDesign() throws URISyntaxException, IOException {
        var resource = getClass().getResource("/stability-design.yaml").toURI();
        var file = new File(resource);
        stabilityDesign = StabilityDesign.parse(file);
    }

    @Test
    void canLoadYAMLFile() {
        Assertions.assertEquals("stability-design", stabilityDesign.designName);
    }

    @Test
    void yamlFilePropertiesAreNotNull() {
        Stream.of(stabilityDesign.designName, stabilityDesign.molecule,
                stabilityDesign.ospreyVersion, stabilityDesign.molecule).forEach(Assertions::assertNotNull);

        Assertions.assertEquals(1, stabilityDesign.molecule.residueModifiers.size());
        Assertions.assertEquals(0.63, stabilityDesign.epsilon, 1e-6);
    }

    @Test
    void yamlFilePdbCanBeParsed() {
        var molecule = PDBIO.read(stabilityDesign.molecule.coordinates);
        assertNotNull(molecule);
    }
}
