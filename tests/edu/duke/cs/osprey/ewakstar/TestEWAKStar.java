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
        double ew = cfp.getParams().getDouble("EW");
        EWAKRatios ewr = new EWAKRatios(cfp);
    }
    
    private ConfigFileParser make2RL0Config() {

        // read config from files
        ConfigFileParser cfp = new ConfigFileParser(new String[]{
            "-c",
            "test/2RL0.kstar/cfgKStar.txt",
            "Dummy command",
            "test/2RL0.kstar/cfgMutSearch.txt",
            "test/2RL0.kstar/cfgSystem.txt"
        });
        cfp.loadData();

        // override file-based config
        // I'm guessing most people have at least two cores, so compute the energy matrix a bit faster
        cfp.getParams().setValue("EmatThreads", "2");

        // this test takes several minutes at the config file's value of e=0.95,
        // but it only takes about two minutes at e=0.99
        cfp.getParams().setValue("epsilon", "0.99");

        return cfp;
    }

    @Test
    public void test3Strands() {
        ConfigFileParser cfp = make2RL0Config();
        cfp.getParams().setValue("MinimizationThreads", "2");
        testEWAKStar(cfp);
    }

}
