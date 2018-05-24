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

package edu.duke.cs.osprey.gmec;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.io.File;
import java.util.Arrays;

import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.gmec.GMECFinder;

/**
 *
 * Testing full computation of GMEC
 * 
 * @author mhall44
 */
public class TestFindGMEC extends TestBase {
    
    @Test
    public void test1CC8(){
        
        // IMPORTANT: need to delete cached files for each test to keep the caches from hiding errors
        for (String name : Arrays.asList("1CC8.EMAT.dat", "1CC8.EPICMAT.dat", "1CC8.TUPEXPEMAT.dat")) {
            new File("examples/1CC8.junit/" + name).delete();
        }
        
        // test once without cached files
        test1CC8Gmec();
        
        // test again with cached files
        test1CC8Gmec();
    }
    
    private void test1CC8Gmec() {
        
        //Here's a 7-residue test using EPIC and LUTE with continuous sidechain flexibility
        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(
            "examples/1CC8.junit/KStar.cfg",
            "examples/1CC8.junit/System.cfg",
            "examples/1CC8.junit/DEE.cfg"
        );
        cfp.loadData();
        
        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
        EnergiedConf gmec = gf.calcGMEC().get(0);
        
        //GMECEnergy should be about -70.617
        //Numerical/fitting error could alter it within ~0.1 kcal/mol.
        //If you are running this test after improving the minimizer maybe you'll
        //find something better, but the energy for this conf (and thus the GMEC)
        //should be at least this good
        assertThat(gmec.getEnergy(), isRelatively(-70.617, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {5, 7, 12, 5, 0, 7, 4}));
        /*
        Example GMEC line in confs.txt:
        0 CONF: 5 7 12 5 0 7 4 RESTYPES: Ile Ser Met Glu Ala Gln Leu ROTS: 5 7 12 5 -1 7 4 Lower bound/enumeration energy: -70.62208670920121 Energy: -70.6172123578607 Best so far: -70.61721235798836 EPIC energy: -70.62518856374005
        */
    }
    
    @Test
    public void testDEEPer(){
        
        // IMPORTANT: need to delete cached files for each test to keep the caches from hiding errors
        for (String name : Arrays.asList("1CC8.EMAT.dat", "1CC8.EPICMAT.dat", "1CC8.TUPEXPEMAT.dat")) {
            new File("examples/1CC8.deeper/" + name).delete();
        }
        
        //Here's a 4-residue test using EPIC, LUTE, and DEEPer together.
        //There are two overlapping backrubs
        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(
            "examples/1CC8.deeper/KStar.cfg",
            "examples/1CC8.deeper/System.cfg",
            "examples/1CC8.deeper/DEE.cfg"
        );
        cfp.loadData();
        
        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
        EnergiedConf gmec = gf.calcGMEC().get(0);
        
        assertThat(gmec.getEnergy(), isRelatively(-49.946, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {5, 7, 36, 5}));
        /*
        Example GMEC line in confs.txt:
        0 CONF: 5 7 36 5 RESTYPES: Ile Ser LEU Glu ROTS: 5 7 3 5 Lower bound/enumeration energy: -49.83992986938508 Energy: -49.94563050249759 Best so far: -49.94563050249759 EPIC energy: -49.927475207707005
        */
    }
}
