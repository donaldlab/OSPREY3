/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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




