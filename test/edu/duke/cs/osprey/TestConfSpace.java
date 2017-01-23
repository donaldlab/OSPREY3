/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey;

import java.util.ArrayList;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;

/**
 *
 * Some tests to make sure the conformational space is set up correctly
 * 
 * @author mhall44
 */
public class TestConfSpace extends TestBase {
    
    
    @BeforeClass
    public static void before() {
            initDefaultEnvironment();
    }
    
    @Test
    public void testConfSpaceGeneration(){
        //Testing setup of a ConfSpace object
        
        ArrayList<String> flexibleRes = new ArrayList<>();
        ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
        
        flexibleRes.add("22");//wt: val (1 dihedral, 3 rot)
        flexibleRes.add("24");//lys (4 dihedrals, 27 rot)
        flexibleRes.add("26");//leu (2 dihedrals, 5 rot)
        flexibleRes.add("27");//thr  (2 dihedrals, 18 rot)
        
        for(int r=0; r<4; r++)
            allowedAAs.add(new ArrayList<String>());
        
        //let's start with just wild-type
        
        ConfSpace cs = new ConfSpace("examples/1CC8/1CC8.ss.pdb", flexibleRes, allowedAAs, true, new ArrayList<>(),
                true, new DEEPerSettings(), new ArrayList<>(), new ArrayList<>(), false, false, null);
        
        //assert some things about the space
        //these are based on the Lovell Rotamer library
        assert cs.confDOFs.size() == 9;//all the dihedrals
        assert cs.mutDOFs.size()==4;//one for each residue
        assert cs.numPos==4;
        assert cs.posFlex.get(1).RCs.size()==27;
        assert cs.posFlex.get(1).res.template.name.equalsIgnoreCase("LYS");
        RC rot1 = cs.posFlex.get(1).RCs.get(1);
        assert rot1.rotNum==1;
        assert rot1.DOFs.size()==4;
        assert Math.round( rot1.DOFmax.get(0) ) == 71;
        assert Math.round( rot1.DOFmin.get(0) ) == 53;
        
        //try adding additional AA types
        allowedAAs.get(0).add("PHE");//2 dihedrals, 4 rotamers.  So now two dihedrals needed for res 22
        allowedAAs.get(1).add("ALA");//no dihedrals or rotamers (but should add one RC)

        cs = new ConfSpace("examples/1CC8/1CC8.ss.pdb", flexibleRes, allowedAAs, true, new ArrayList<>(),
                false, new DEEPerSettings(), new ArrayList<>(), new ArrayList<>(), false, false, null);
        
        assert cs.confDOFs.size()==10;
        assert cs.mutDOFs.size()==4;
        assert cs.posFlex.get(0).RCs.size()==7;
        assert cs.posFlex.get(1).RCs.size()==28;
        
        RC rot0 = cs.posFlex.get(2).RCs.get(0);
        assert rot0.rotNum==0;
        assert rot0.DOFs.size()==2;
        assert Math.round( rot0.DOFmax.get(1) ) == 80;
        assert Math.round( rot0.DOFmin.get(1) ) == 80;
        
        //if we get here, test passed
        System.out.println("CONFSPACE GENERATION TEST PASSED");
    }
}
