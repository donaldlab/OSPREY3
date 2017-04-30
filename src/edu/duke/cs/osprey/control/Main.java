/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.energy.LigandResEnergies;
import java.util.HashMap;
import java.util.Map;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import edu.duke.cs.osprey.tests.UnitTestSuite;
import edu.duke.cs.osprey.tools.Stopwatch;

/**
 *
 * @author mhall44
 * Parse arguments and call high-level functions like DEE/A* and K*
   These will each be controlled by dedicated classes, unlike in the old KSParser
   to keep this class more concise, and to separate these functions for modularity purposes
 */

public class Main {

	public static Map<String, Runnable> commands;
        
        private static final String usageString = "Command expects arguments "
                + "(e.g. -c KStar.cfg {findGMEC|calcKStar} System.cfg DEE.cfg";

	public static void main(String[] args){
		//args expected to be "-c KStar.cfg command config_file_1.cfg ..."

		debuggingCommands(args);

		String command = "";
		try{
                    command = args[2];
		}
		catch(Exception e){
			System.out.println(usageString);
			System.exit(1);
		}



		Stopwatch stopwatch = new Stopwatch().start();

		ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files

                EnvironmentVars.openSpecialWarningLogs(cfp);

		//load data files
		cfp.loadData();



		//DEBUG!!
		initCommands(args, cfp);

		if(commands.containsKey(command))
			commands.get(command).run();
		else
			throw new RuntimeException("ERROR: OSPREY command unrecognized: "+command);
                
                EnvironmentVars.closeSpecialWarningLogs();
                
		System.out.println("Total OSPREY execution time: " + stopwatch.getTime(2));
		System.out.println("OSPREY finished");
	}

	private static void initCommands(String[] args, ConfigFileParser cfp) {
		
		// set degree of thread parallelism
		// NOTE: if we're going to use the config files here, don't override its defaults
		ThreadParallelism.setNumThreads(cfp.params.getInt("NumThreads"));
		MultiTermEnergyFunction.setNumThreads(ThreadParallelism.getNumThreads());

                CCDMinimizer.EConvTol = cfp.params.getDouble("CCDEConvTol");
                CCDMinimizer.numIter = cfp.params.getInt("CCDNumIter");
                EnvironmentVars.alwaysIdealizeSidechainsAfterMutation = cfp.params.getBool("ALWAYSIDEALIZESIDECHAINSAFTERMUTATION");
                
		// TODO Auto-generated method stub
		commands = new HashMap<String, Runnable>();

		commands.put("findGMEC", new Runnable() {
			@Override
			public void run() {
				GMECFinder gf = new GMECFinder();
				gf.init(cfp);
				gf.calcGMEC();
			}
		});
		
		commands.put("findSequences", new Runnable() {
			@Override
			public void run() {
				GMECFinder gf = new GMECFinder();
				gf.init(cfp);
				gf.calcSequences();
			}
		});
		
		commands.put("scoreLowestConfs", new Runnable() {
			@Override
			public void run() {
				BigForcefieldEnergy.ParamInfo.printWarnings = false;
				ForcefieldParams.printWarnings = false;
				
				GMECFinder gf = new GMECFinder();
				gf.init(cfp);
				gf.setLogConfsToConsole(false);
				gf.scoreLowestConfs();
			}
		});

		commands.put("calcKStar", new Runnable() {
			@Override
			public void run() {
				// kstar subclasses configfileparser, so re-load
				KSConfigFileParser ksCfp = new KSConfigFileParser(args);
				ksCfp.loadData();
				
				KStarCalculator ksc = new KStarCalculator(ksCfp);
				ksc.calcKStarScores();
			}
		});

		commands.put("RunTests", new Runnable() {
			@Override
			public void run() {
				UnitTestSuite.runAllTests();
			}
		});

		commands.put("doCOMETS", new Runnable() {
			@Override
			public void run() {
				COMETSDoer cd = new COMETSDoer(args);
				cd.calcBestSequences();
			}
		});
		
		commands.put("doMSKStar", new Runnable() {
			@Override
			public void run() {
				MSKStarDoer msksd = new MSKStarDoer(args);
				msksd.calcBestSequences();
			}
		});

		commands.put("calcLigResE", new Runnable() {
			@Override
			public void run() {
				LigandResEnergies lre = new LigandResEnergies(cfp.getParams());
				lre.printEnergies();
			}
		});

		commands.put("calcEnergy", new Runnable() {
			@Override
			public void run() {
				new EnergyCalculator().run(cfp);
			}
		});

		commands.put("ConfInfo", new Runnable() {
			@Override
			public void run() {
				ConfInfo ci = new ConfInfo(cfp);
				ci.outputConfInfo();
			}
		});
                
                commands.put("findSeqGMECs", new Runnable() {
			@Override
			public void run() {
				SeqGMECFinder sgf = new SeqGMECFinder(args, cfp.getParams().getValue("MutFile"));
                                sgf.calcAllSeqGMECs();
			}
		});
	}

	// TODO: Move these into a test file, and just call it from the test.
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
