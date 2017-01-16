/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.GMECFinder;

/**
 *
 * Testing full computation of GMEC (precomputed energy matrices will be used,
 * but delete them if you want to test that part too)
 * 
 * @author mhall44
 */
public class TestFindGMEC extends TestBase {
    
    @Test
    public void test1CC8(){
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
