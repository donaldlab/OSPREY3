/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

/**
 *
 * @author mhall44
 */
public class Main {
    //Parse arguments and call high-level functions
    //like DEE/A* and K*
    //These will each be controlled by dedicated classes, unlike in the old KSParser
    //to keep this class more concise, and to separate these functions for modularity purposes
    
    /*TODO (as of 3/26/15):
    1. Figure out template structure & assignment
    2. Figure out energy/pruning matrix structure including higher-order
    3. Pass through everything and try to handle incomplete stuff, except what is to be deferred til after gotten working
    4. List stuff for deferral, comment out
    5. Get working, debug/test, and show everybody
    6. Make changes based on comments and distribute to people who want to help
    7. Implement additional features and concurrently start building up tuple expander (in very modular way
    * so don't have to change rest of code), etc.  
    (Tuple expander can go directly for two goals: huge designs, Poisson-Boltzmann.  Entropy, BB after that).  
    */
    
    
    public static void main(String[] args){
        //args expected to be "-c KStar.cfg command config_file_1.cfg ..."
        String command = args[2];
        long startTime = System.currentTimeMillis();
        
        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
        
        //load data files
        cfp.loadData();
        
        if(command.equalsIgnoreCase("findGMEC")){
            //I recommend that we change the command names a little to be more intuitive, e.g. 
            //"findGMEC" instead of doDEE
            GMECFinder gf = new GMECFinder(cfp);
            gf.calcGMEC();//this can be the n globally minimum-energy conformations for n>1, or just the 1 
            //These functions will handle their own output
        }
        else if(command.equalsIgnoreCase("calcKStar")){
            throw new UnsupportedOperationException("ERROR: Still need to implement K*");
            //KStarCalculator ksc = new KStarCalculator(params);
            //ksc.calcKStarScores();
        }
        //etc.
        else
            throw new RuntimeException("ERROR: OSPREY command unrecognized: "+command);
        
        long endTime = System.currentTimeMillis();
        System.out.println("Total OSPREY execution time: " + ((endTime-startTime)/60000) + " minutes.");
        System.out.println("OSPREY finished");
    }
    
    
    
}
