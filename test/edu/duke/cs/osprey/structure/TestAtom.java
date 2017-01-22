package edu.duke.cs.osprey.structure;

import static org.junit.Assert.*;
import static org.hamcrest.Matchers.*;

import org.junit.Test;

public class TestAtom {
	
	@Test
	public void nameToElement() {
		assertElement("C", "C");
		assertElement("H", "H");
		assertElement("Zn", "Zn");
		assertElement("CA", "C"); // not calcium, the A stands for alpha in PDB names
		assertElement("CB", "C"); // C beta
		assertElement("CG", "C"); // C gamma
		assertElement("CD", "C"); // C delta
		assertElement("CE", "C"); // C epsilon
		assertElement("Ca", "Ca"); // calcium
		assertElement("Cd", "Cd"); // cadmium
		assertElement("Ce", "Ce"); // cesium
		assertElement("NA", "N"); // N alpha
		assertElement("Na", "Na"); // sodium
		
		// can't recognize zinc without proper capitalization
		// if we assume a case for atom names, we'll break detection for e.g. CA (alpha carbon) vs Ca (calcium)
		assertElement("ZN", "DU");
	}
	
	private void assertElement(String atomName, String elem) {
		assertThat(new Atom(atomName).elementType, is(elem));
	}
}
