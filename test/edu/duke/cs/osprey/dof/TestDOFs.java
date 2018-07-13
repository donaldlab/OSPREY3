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
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
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

package edu.duke.cs.osprey.dof;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

/**
 *
 * @author mhall44
 */
public class TestDOFs extends TestBase {
    
    @Test
    public void testMutation() {
        
        Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
        Residue res = strand.mol.residues.get(37); // Ser 39 originally
        
        res.pucker = new ProlinePucker(strand.templateLib, res);
        ResidueTypeDOF mutDOF = new ResidueTypeDOF(strand.templateLib, res);
        
        mutDOF.mutateTo("ALA");
        assertThat(res.template.name, is("ALA"));
        
        mutDOF.mutateTo("ILE");
        assertThat(res.template.name, is("ILE"));
        
        mutDOF.mutateTo("VAL");
        assertThat(res.template.name, is("VAL"));
        
        mutDOF.mutateTo("PRO");
        assertThat(res.template.name, is("PRO"));
        
        mutDOF.mutateTo("GLY");
        assertThat(res.template.name, is("GLY"));
        
        mutDOF.mutateTo("ARG");
        assertThat(res.template.name, is("ARG"));
    }
    
    @Test
    public void testDihedral(){
        
        Molecule m = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build().mol;
        Residue res = m.residues.get(37);
        
        FreeDihedral chi1 = new FreeDihedral(res,0);//Ser 39
        FreeDihedral chi2 = new FreeDihedral(res,1);
        
        chi1.apply(45);
        chi2.apply(-121);
        
        //measure dihedrals.  Start by collecting coordinates
        double N[] = res.getCoordsByAtomName("N");
        double CA[] = res.getCoordsByAtomName("CA");
        double CB[] = res.getCoordsByAtomName("CB");
        double OG[] = res.getCoordsByAtomName("OG");
        double HG[] = res.getCoordsByAtomName("HG");
        
        assertThat(Protractor.measureDihedral(new double[][] {N,CA,CB,OG}), isRelatively(45));
        assertThat(Protractor.measureDihedral(new double[][] {CA,CB,OG,HG}), isRelatively(-121));
    }
}
