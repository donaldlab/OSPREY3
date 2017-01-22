package edu.duke.cs.osprey;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.List;
import java.util.Map;

import org.junit.Test;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.restypes.PositionSpecificRotamerLibrary;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;

public class TestGenerateRotamerLibrary {
    
    @Test
    public void testLoadPDBRotamers()
    throws Exception {
    	
        Strand strand = Strand.builder(PDBIO.readFile("examples/4NPD/4NPD.pdb"))
            .setResidues(1, 58)
            .build();
        
        PositionSpecificRotamerLibrary lib = new PositionSpecificRotamerLibrary(strand);
        
        // check the templates and rotamers
        for (int i=0; i<strand.mol.residues.size(); i++) {
        	
        	Residue res = strand.mol.residues.get(i);
        	Map<String,List<ResidueTemplate>> templatesByAA = lib.getTemplatesForDesignIndex(i);
        	
        	assertThat(templatesByAA.keySet(), contains(res.template.name));
        	assertThat(templatesByAA.get(res.template.name).size(), is(1));
        	assertThat(templatesByAA.get(res.template.name).get(0).getNumRotamers(0, 0), is(1 + strand.mol.getAlternates(i).size()));
        }
    }
}
