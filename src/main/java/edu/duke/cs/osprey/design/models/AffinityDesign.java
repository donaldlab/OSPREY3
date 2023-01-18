package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.*;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class AffinityDesign {
    @JsonProperty("osprey_version")
    public String ospreyVersion = "3.1.0";

    @JsonProperty("design_name")
    public String designName = "default-design";

    @JsonProperty("protein")
    public MoleculeDto protein;

    @JsonProperty("ligand")
    public MoleculeDto ligand;

    @JsonProperty("scan")
    public ScanDto scanSettings;

    @JsonProperty("sequence_filters")
    public List<String> sequenceFilters;

    private static final YAMLFactory factory = new YAMLFactory()
            .enable(YAMLGenerator.Feature.MINIMIZE_QUOTES)
            .enable(YAMLGenerator.Feature.LITERAL_BLOCK_STYLE);

    public static AffinityDesign parse(File file) throws IOException {
        var mapper = new ObjectMapper(factory);
        mapper.configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, false);
        var stream = new FileInputStream(file);
        return mapper.readValue(stream, AffinityDesign.class);
    }

    public void write(OutputStream stream) throws IOException {
        var mapper = new ObjectMapper(factory);
        mapper.setSerializationInclusion(JsonInclude.Include.NON_NULL);
        mapper.writeValue(stream, this);
    }

    public void write(Path dest) throws IOException {
        var stream = new FileOutputStream(dest.toFile());
        write(stream);
    }

    public AffinityDesign copy() {
        var om = new ObjectMapper();
        try {
            return om.readValue(om.writeValueAsString(this), this.getClass());
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    public List<String> validate() {
        var errors = new ArrayList<String>();

        if (this.ligand == null || this.protein == null) {
            errors.add("The this must have both a protein and ligand specified");
            return errors;
        }

        if (this.protein.coordinates.isEmpty()) {
            errors.add("You must specify PDB coordinates for the protein");
        }

        if (this.ligand.coordinates.isEmpty()) {
            errors.add("You must specify PDB coordinates for the ligand");
        }

        return errors;
    }

    public Molecule makeProteinMolecule() {
        return PDBIO.read(protein.coordinates);
    }

    public Molecule makeLigandMolecule() {
        return PDBIO.read(ligand.coordinates);
    }
}
