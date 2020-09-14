package edu.duke.cs.osprey.design.models;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
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

    @JsonProperty("epsilon")
    public double epsilon = 0.63;

    @JsonProperty("scan")
    public ScanDto scanSettings;

    public static AffinityDesign parse(File file) throws IOException {
        var mapper = new ObjectMapper(new YAMLFactory());
        var stream = new FileInputStream(file);
        return mapper.readValue(stream, AffinityDesign.class);
    }

    public void write(Path dest) throws IOException {
        var mapper = new ObjectMapper(new YAMLFactory());
        mapper.setSerializationInclusion(JsonInclude.Include.NON_NULL);
        var stream = new FileOutputStream(dest.toFile());
        mapper.writeValue(stream, this);
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

        if (this.epsilon <= 0 || this.epsilon >= 1) {
            errors.add(String.format("Epsilon must be between 0 and 1 (exclusive). Was %f", this.epsilon));
        }

        if (this.ligand == null || this.protein == null) {
            errors.add("The this must have both a protein and ligand specified");
        }

        if (this.protein.coordinates.isEmpty()) {
            errors.add("You must specify PDB coordinates for the protein");
        }

        if (this.ligand.coordinates.isEmpty()) {
            errors.add("You must specify PDB coordinates for the ligand");
        }

        if (this.protein.residueModifiers.size() + this.ligand.residueModifiers.size() == 0) {
            errors.add("There must be at least one flexible/mutable residue in the this");
        }

        return errors;
    }
}
