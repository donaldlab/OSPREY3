/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests.variationalkstar;

import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.SublinearKStarDoer;
import edu.duke.cs.osprey.tools.Stopwatch;
import java.io.File;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author hmn5
 */
public class SublinearKStarTests {

    public static void main(String[] args)
            throws Exception {

        String path = new File("").getAbsolutePath();
        if (!path.endsWith("VariationalKStar/SublinearKStar/4HEM")) {
            throw new Error("This test was designed to be run in testVariationalKStar/SublinearKStar/4HEM folder\n\tcwd: " + path);
        }

        //load configurations
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();
        
        PartFuncTree.verbose = false;
        SublinearKStarDoer skd = new SublinearKStarDoer(cfp);
        if (args[2].equalsIgnoreCase("DOEXHAUSTIVE")) {
            Stopwatch.start();
            skd.doSublinearKStar(true);
            Stopwatch.stop();
            System.out.println("Exhaustive took: "+Stopwatch.getTime(TimeUnit.MILLISECONDS));
        } else {
            Stopwatch.start();
            skd.doSublinearKStar(false);
            Stopwatch.stop();
            System.out.println("Sublinear KStar took: "+Stopwatch.getTime(TimeUnit.MILLISECONDS));
        }
    }
    
}
