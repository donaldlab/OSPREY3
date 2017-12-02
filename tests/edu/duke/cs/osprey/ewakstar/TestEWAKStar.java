/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.ewakstar.EWAKRatios;

import org.junit.Test;

/**
*
* @author Anna Lowegard(anna.lowegard@duke.edu)
* Adegoke Ojewole (ao68@duke.edu)
*/
public class TestEWAKStar {
	
	public void testEWAKStar(ConfigFileParser cfp) {
        EWAKRatios ewr = new EWAKRatios(cfp);
        ewr.run();
    }
    
    private ConfigFileParser make2RL0Config() {

        // read config from files
        ConfigFileParser cfp = new ConfigFileParser(new String[]{
            "-c",
            "test/2RL0.ewakstar/cfgKStar.txt",
            "Dummy command",
            "test/2RL0.ewakstar/cfgMutSearch.txt",
            "test/2RL0.ewakstar/cfgSystem.txt"
        });
        cfp.loadData();    
        return cfp;
    }

    @Test
    public void test3Strands() {
        ConfigFileParser cfp = make2RL0Config();
        // override file-based config
        cfp.getParams().setValue("MinimizationThreads", "2");
        cfp.getParams().setValue("EmatThreads", "2");
        cfp.getParams().setValue("IVAL", "5.0");
        cfp.getParams().setValue("EW", "5.0");
        testEWAKStar(cfp);
    }
    
    @Test
    public void test3StrandsAndMutFile() {
        ConfigFileParser cfp = make2RL0Config();
        // override file-based config
        cfp.getParams().setValue("MinimizationThreads", "2");
        cfp.getParams().setValue("EmatThreads", "2");
        cfp.getParams().setValue("IVAL", "0.0");
        cfp.getParams().setValue("EW", "0.0");
        cfp.getParams().setValue("MUTFILE", "test/2RL0.ewakstar/mutFile.txt");
        testEWAKStar(cfp);
    }

}
