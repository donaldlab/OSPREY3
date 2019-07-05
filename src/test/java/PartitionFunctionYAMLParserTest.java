import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import edu.duke.cs.osprey.design.models.StabilityDesign;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Assertions;

import java.io.IOException;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertNotNull;

class PartitionFunctionYAMLParserTest {

    private StabilityDesign stabilityDesign;

    @BeforeEach
    void getTestDesign() {
        var mapper = new ObjectMapper(new YAMLFactory());
        var stream = getClass().getResourceAsStream("test-design.yaml");

        try {
            stabilityDesign = mapper.readValue(stream, StabilityDesign.class);
        } catch (IOException e) {
            Assertions.fail(e);
        }
    }

    @Test
    void canLoadYAMLFile() {
        Assertions.assertEquals("Test design", stabilityDesign.designName);
    }

    @Test
    void yamlFilePropertiesAreNotNull() {
        Stream.of(stabilityDesign.designName, stabilityDesign.molecule,
                stabilityDesign.ospreyVersion, stabilityDesign.residueModifiers).forEach(Assertions::assertNotNull);

        Assertions.assertEquals(8, stabilityDesign.residueModifiers.size());
        Assertions.assertEquals(0.63, stabilityDesign.epsilon, 1e-6);
    }

    @Test
    void yamlFilePdbCanBeParsed() {
        var molecule = PDBIO.read(stabilityDesign.molecule);
        assertNotNull(molecule);
    }
}
