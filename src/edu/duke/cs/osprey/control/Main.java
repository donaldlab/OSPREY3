/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.tests.UnitTestSuite;

/**
 *
 * @author mhall44
 * Parse arguments and call high-level functions like DEE/A* and K*
   These will each be controlled by dedicated classes, unlike in the old KSParser
   to keep this class more concise, and to separate these functions for modularity purposes
 */

public class Main {

    
    
    public static void main(String[] args){
        //args expected to be "-c KStar.cfg command config_file_1.cfg ..."
        
        debuggingCommands(args);
        
        String command = "";
        try{
        	command = args[2];
        }
        catch(Exception e){
        	System.out.println("Command expects arguments (e.g. -c KStar.cfg {findGMEC|fcalcKStar} System.cfg DEE.cfg");
        	System.exit(1);
        }
        
        
        
        long startTime = System.currentTimeMillis();
        
        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
        
        //load data filescloneclone
        cfp.loadData();        
        
        
        //DEBUG!!
        // set number of threads for energy function evaluation
        MultiTermEnergyFunction.setNumThreads( cfp.params.getInt("eEvalThreads", 1) );
        if( MultiTermEnergyFunction.getNumThreads() > 1 ) {
                System.setProperty( "java.util.concurrent.ForkJoinPool.common.parallelism", 
                                String.valueOf(MultiTermEnergyFunction.getNumThreads()) );
        }
        
        
        
        if(command.equalsIgnoreCase("findGMEC")){
            //I recommend that we change the command names a little to be more intuitive, e.g. 
            //"findGMEC" instead of doDEE
            GMECFinderSuper gf = new GMECFinderSuper(cfp);
            gf.calcGMEC();//this can be the n globally minimum-energy conformations for n>1, or just the 1 
            //These functions will handle their own output
        }
        else if(command.equalsIgnoreCase("calcKStar")){
            throw new UnsupportedOperationException("ERROR: Still need to implement K*");
            //KStarCalculator ksc = new KStarCalculator(params);
            //ksc.calcKStarScores();
        }
        else if(command.equalsIgnoreCase("RunTests")){
            UnitTestSuite.runAllTests();
        }
        else if(command.equalsIgnoreCase("doCOMETS")){
            COMETSDoer cd = new COMETSDoer(args);
            cd.calcBestSequences();
        }
        else if(command.equalsIgnoreCase("KaDEETest")){
            KaDEEFinder kdf = new KaDEEFinder(cfp);
            kdf.doKaDEE();
        }
        //etc.
        else
            throw new RuntimeException("ERROR: OSPREY command unrecognized: "+command);
        
        long endTime = System.currentTimeMillis();
        System.out.println("Total OSPREY execution time: " + ((endTime-startTime)/60000) + " minutes.");
        System.out.println("OSPREY finished");
    }
    
    
    private static void debuggingCommands(String[] args){
        //anything we want to try as an alternate main function, for debugging purposes
        //likely will want to exit after doing this (specify here)
        //for normal operation, leave no uncommented code in this function
        
    }
    
    
}