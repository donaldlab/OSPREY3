package edu.duke.cs.osprey.design;

import edu.duke.cs.osprey.design.models.AffinityDesign;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.List;

class AffinityDesignYAMLParserTest {

    private AffinityDesign affinityDesign;


    @BeforeEach
    void getTestDesign() throws IOException, URISyntaxException {
        affinityDesign = getDesignInResources("/affinity-design.yaml");
    }

    @Test
    void canLoadYAMLFile() {
        Assertions.assertEquals("affinity-design", affinityDesign.designName);
    }

    @Test
    void hasProteinAndLigand() {
        Assertions.assertNotNull(affinityDesign.protein);
        Assertions.assertNotNull(affinityDesign.ligand);
        Assertions.assertNotNull(affinityDesign.protein.coordinates);
        Assertions.assertNotNull(affinityDesign.ligand.coordinates);
    }

    @Test
    void passesValidation() {
        Assertions.assertIterableEquals(List.of(), affinityDesign.validate());
    }

    @Test
    void failsValidation() throws Exception {
        var design = getDesignInResources("/bad-affinity-design.yaml");
        Assertions.assertNotEquals(0, design.validate().size());
    }

    private AffinityDesign getDesignInResources(String named) throws IOException, URISyntaxException {
        var resource = getClass().getResource(named).toURI();
        var file = new File(resource);
        return AffinityDesign.parse(file);
    }
}
