/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests.variationalkstar;

import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.partitionfunctionbounds.SequenceFreeEnergy;
import edu.duke.cs.osprey.partitionfunctionbounds.TRBP;
import edu.duke.cs.osprey.tools.Stopwatch;
import java.io.File;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author hmn5
 */
public class SublinearFreeEnergyTests {

    public static void main(String[] args)
            throws Exception {

        String path = new File("").getAbsolutePath();
        if (!path.endsWith("VariationalKStar/SublinearFreeEnergy/4HEM")) {
            throw new Error("This test was designed to be run in testVariationalKStar/SublinearFreeEnergy/4HEM folder\n\tcwd: " + path);
        }

        //load configurations
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();

        SequenceFreeEnergy sf = new SequenceFreeEnergy(cfp);

        if (args[2].equalsIgnoreCase("DOEXHAUSTIVE")) {
            sf.exhaustiveSearch();
        }

        Stopwatch.start();
        int[] bestSequence = sf.nextConf();
        String[] seqList = sf.bestSequence;
        Stopwatch.stop();
        String correctSequence = "ARG ARG ARG ARG";
        
        StringBuilder sb = new StringBuilder();
        for (String aa : seqList) {
            sb.append(aa);
            sb.append(" ");
        }
        String sequence = sb.toString();

        System.out.println("Finished in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));
        System.out.println();
        System.out.println("Sublinear Sequence:  " + sequence);
        System.out.println("Correct Sequence  :  " + correctSequence);
    }
}
