/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.util.HashMap;
import java.util.Map;
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
        MultiTermEnergyFunction.setNumThreads( cfp.params.getInt("eEvalThreads") );
        if( MultiTermEnergyFunction.getNumThreads() > 1 ) {
                System.setProperty( "java.util.concurrent.ForkJoinPool.common.parallelism", 
                                String.valueOf(MultiTermEnergyFunction.getNumThreads()) );
        }
        
        
        
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
        else if(command.equalsIgnoreCase("RunTests")){
            UnitTestSuite.runAllTests();
        }
        else if(command.equalsIgnoreCase("doCOMETS")){
            COMETSDoer cd = new COMETSDoer(args);
            cd.calcBestSequences();
//            cd.exhaustiveMultistateSearch();
        }
        else if(command.equalsIgnoreCase("doCOMETSSuper")){
            COMETSDoerSuper csd = new COMETSDoerSuper(cfp);
            csd.doKaDEE();
        }
        else if (command.equalsIgnoreCase("doKaDEE")){
            KaDEEDoer kd = new KaDEEDoer(cfp);
            kd.doKaDEE();
        }
        else if (command.equalsIgnoreCase("doGumbel")){
            GumbelDoer gd = new GumbelDoer(cfp);
//            KaDEEDoer2 kd = new KaDEEDoer2(cfp);
//            kd.doKaDEE();
        }
        //etc.
        else
            throw new RuntimeException("ERROR: OSPREY command unrecognized: "+command);
        
        long endTime = System.currentTimeMillis();
        System.out.println("Total OSPREY execution time: " + ((endTime-startTime)/60000) + " minutes.");
        System.out.println("OSPREY finished");
    }
    
    
    private static void debuggingCommands(String[] args){
        
		//MolecEObjFunction mof = (MolecEObjFunction)ObjectIO.readObject("OBJFCN1442697734046.dat", true);
		/*MolecEObjFunction mof = (MolecEObjFunction)ObjectIO.readObject("OBJFCN1442697735769.dat", true);

        CCDMinimizer minim = new CCDMinimizer(mof, false);
        DoubleMatrix1D bestVals = minim.minimize();
        double E = mof.getValue(bestVals);

        DoubleMatrix1D boxBottom = bestVals.copy();
        DoubleMatrix1D boxTop = bestVals.copy();
        for(int q=0; q<bestVals.size(); q++){
            boxTop.set(q, Math.min(mof.getConstraints()[1].get(q), bestVals.get(q)+1));
            boxBottom.set(q, Math.max(mof.getConstraints()[0].get(q), bestVals.get(q)-1));
        }

        for(int a=0; a<1000000; a++){

            DoubleMatrix1D x2 = bestVals.like();
            for(int q=0; q<bestVals.size(); q++)
                x2.set( q, boxBottom.get(q)+Math.random()*(boxTop.get(q)-boxBottom.get(q)) );

            double E2 = mof.getValue(x2);
            if(E2 < E-0.1){
                System.out.println("gg");
                DoubleMatrix1D dx = x2.copy();
                dx.assign( bestVals, Functions.minus );

                for(double t=1; true; t*=1.5){
                    dx.assign(Functions.mult(t));
                    x2.assign(dx);
                    x2.assign(bestVals, Functions.plus);

                    boolean oor = false;
                    for(int q=0; q<x2.size(); q++){
                        if( x2.get(q) > mof.getConstraints()[1].get(q) )
                            oor = true;
                        else if( x2.get(q) < mof.getConstraints()[0].get(q) )
                            oor = true;
                    }

                    if(oor)
                        break;

                    E2= mof.getValue(x2);
                    int aaa = 1;
                }
            }
        }

        System.exit(0);
		 */
		//anything we want to try as an alternate main function, for debugging purposes
		//likely will want to exit after doing this (specify here)
		//for normal operation, leave no uncommented code in this function

	}


}
