package edu.duke.cs.osprey.triplesBounds;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;


public class TestTriplesEnergy {

	@Test
	public void indexing() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A39").setLibraryRotamers("VAL").setContinuous();
		strand.flexibility.get("A40").setLibraryRotamers("ILE").setContinuous();
		strand.flexibility.get("A41").setLibraryRotamers("LEU").setContinuous();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

		TriplesEnergy energies = new TriplesEnergy(confSpace);
		energies.fill(0.0);

		for (int pos1=0; pos1<confSpace.positions.size(); pos1++) {
			int numRCs1 = confSpace.positions.get(pos1).resConfs.size();
			for (int pos2=0; pos2<pos1; pos2++) {
				int numRCs2 = confSpace.positions.get(pos2).resConfs.size();
				for (int pos3=0; pos3<pos2; pos3++) {
					int numRCs3 = confSpace.positions.get(pos3).resConfs.size();

					for (int rc1=0; rc1<numRCs1; rc1++) {
						for (int rc2=0; rc2<numRCs2; rc2++) {
							for (int rc3=0; rc3<numRCs3; rc3++) {

								assertThat(energies.get(pos1, rc1, pos2, rc2, pos3, rc3), is(0.0));
								energies.set(pos1, rc1, pos2, rc2, pos3, rc3, 1.0);
								assertThat(energies.get(pos1, rc1, pos2, rc2, pos3, rc3), is(1.0));
							}
						}
					}
				}
			}
		}
	}
}
