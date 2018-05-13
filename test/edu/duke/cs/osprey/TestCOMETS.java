/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/


package edu.duke.cs.osprey;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.COMETSDoer;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Test;

/**
 *
 * Testing full COMETS computation
 * 
 * @author mhall44
 */
public class TestCOMETS {
    
    @Test
    public void testCOMETS() {
    	
    	// IMPORTANT: delete the cached energy matrices before running this test!
    	// cached energy matrices can hide errors involving config or emat calculation
    	for (String path : Arrays.asList("3K75.b.EMAT.dat", "3K75.ub.EMAT.dat", "3LQC.b.EMAT.dat", "3LQC.ub.EMAT.dat")) {
    		new File("examples/comets.junit/" + path).delete();
    	}
    	
        //Here's a specificity design with discrete flexibility, 5 mutable residues
        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(
            "examples/comets.junit/KStar.cfg",
            "examples/comets.junit/multistate.spec0.cfg"
        );
        cfp.loadData();
        
        COMETSDoer cd = new COMETSDoer(cfp);
        ArrayList<String> bestSequences = cd.calcBestSequences();
        
        
        String[] actualBestSequences = new String[] { "TYR_GLU_ILE_LEU_GLN_", 
            "TYR_GLU_MET_MET_GLN_", "TYR_GLU_PHE_LEU_GLN_",
            "TYR_GLU_VAL_PHE_GLN_", "TYR_GLU_VAL_LEU_GLN_" };
        
        //OK let's make sure the returned 5 best sequences match the actual ones
        //we can expect an exact match here because the search is discrete
        assertThat(bestSequences.size(), is(5));
        for(int seqNum=0; seqNum<5; seqNum++) {
            assertThat(bestSequences.get(seqNum), is(actualBestSequences[seqNum]));
        }
    }
}

