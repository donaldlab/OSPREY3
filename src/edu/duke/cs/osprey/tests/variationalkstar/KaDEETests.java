/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests.variationalkstar;

import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.KaDEEDoer;
import java.io.File;

/**
 *
 * @author hmn5
 */
public class KaDEETests {

    static String[] dirNums = {"01", "02", "03", "04", "05", "06", "07", "08", "09"};
//        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
//        "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"};

    public static void main(String[] args) {
        String path = new File("").getAbsolutePath();
        if (!path.endsWith("VariationalKStar/KaDEE")) {
            throw new Error("This test was designed to be run in test/VariationalKStar/KaDEE folder\n\tcwd: " + path);
        }
        String[] pathsToDir = {"LargeTest/4LAJ/", "LargeTest/3GXU/", "LargeTest/4HEM/"};
//         String[] pathsToDir = {"LargeTest/4HEM/"};
        for (String pathToDir : pathsToDir) {
            for (String dir : dirNums) {
                String runPath = pathToDir + dir + "_Run/";
                ConfigFileParser cfp = setupRun(runPath);
                KaDEEDoer kd = new KaDEEDoer(cfp);
                kd.doKaDEE();
            }
        }
    }

    private static ConfigFileParser setupRun(String run) {
        String dee = run + "DEE.cfg";
        String sys = run + "System.cfg";
        String kstar = run + "KStar.cfg";

        String[] newArgs = {"-c", kstar, " ", dee, sys};
        ConfigFileParser cfp = new ConfigFileParser(newArgs);
        String pdbFile = cfp.params.getValue("PDBNAME");
        String pathToPDB = run + pdbFile;
        cfp.params.setValue("PDBNAME", pathToPDB);
        String name = cfp.params.getValue("RUNNAME");
        String pathName = run + name;
        cfp.params.setValue("RUNNAME", pathName);
        cfp.params.setValue("DOSOLVATIONE", "FALSE");
//        cfp.params.setValue("ADDWTROTS", "FALSE");
        cfp.loadData();

        return cfp;
    }
}
