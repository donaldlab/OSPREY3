/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.gmec;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfPrinter;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.partcr.PartCRConfPruner;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Progress;
import edu.duke.cs.osprey.tools.Stopwatch;

/**
 *
 * @author mhall44
 */
public class GMECFinder {
	
    public static interface ConfPruner {
        void prune(List<ScoredConf> confs, GMECConfEnergyCalculator ecalc);
    }
    
    //Many of the parameters that control the optimization process (like what DEE algs to use,
    //whether to use MPLP, etc. can be stored in this class
    //I propose that they be grouped into classes, similar to "EPICSettings," that handle settings for
    //particular aspects of functionality (energy function parameters, DEE/A* alg settings, etc.)
    //But we can also just store things as fields if people prefer
    
    //KStarCalculator will be set up similarly to this class
        
    private SearchProblem searchSpace;
    private PruningControl pruningControl;
    private ConfSearchFactory confSearchFactory;
    private GMECConfEnergyCalculator.Async ecalc;
    private ConfPruner confPruner;
    
    private double Ew;//energy window for enumerating conformations: 0 for just GMEC
    private double I0=0;//initial value of iMinDEE pruning interval
    private boolean doIMinDEE;//do iMinDEE
    
    
    private boolean useContFlex;
    //boolean enumByPairwiseLowerBound;//are we using a search method that 
    //enumerates by pairwise lower bound (minDEE-style)
    //or by (possibly approximated) true energy (rigid, EPIC, or tuple expander)?
    
    
    private boolean outputGMECStruct;//write the GMEC structure to a PDB file
    
    boolean useEPIC = false;
    public boolean useTupExp = false;
    
    private boolean checkApproxE = true;//Use the actual energy function to evaluate
    //each enumerated conformation rather than just relying on the EPIC or tup-exp approximation
    
    public boolean EFullConfOnly = false;//energy function only can be evaluated for full conf
    
    private String confFileName;//file to which we write conformations
    
    private double stericThresh;
    private boolean logConfsToConsole;
    
    private double lowestBound;

    private boolean usePLUG;//Use PLUG, including multi-term pruning
    static boolean PLUGPruneTriples = true;//Seems fast enough (cubic time--no competitors)
    private boolean useDEE = true;//Sometimes may not want to prune with DEE, esp if using PLUG
    //(note: will skip DEE regardless if not pairwise E-fcn)


    public GMECFinder() {
        
        // Arguably, the stuff below is the initialization, which should be its own function. The constructor may eventually
        // do other more interesting things and reading in member variables isn't that interesting. -JJ
        // Ew = cfp.prams.getDouble("Ew");
        // etc...
    
        // Jeff: I think is a good place to set default values
        logConfsToConsole = true;
        confPruner = null;
        ecalc = null;
        
        // TODO: in the future, we want all config options to act like the ones set in this constructor
        // The goal is: construct a GMECFinder instance without huge arguments lists and it will use defaults automatically
        // then, when you want to deviate from defaults, you call e.g. :  gmecFinder.setThisOption(toThisVal)
        // this will make the eventual python scripting interface much much nicer!
        // options that can't fundamentally be defaulted are actually inputs, and they should generally get passed
        // in the method args in the method that actually needs them
    }
    
    public void init(SearchProblem search, PruningControl pruningControl, ConfSearchFactory confSearchFactory, GMECConfEnergyCalculator.Async ecalc,
        double Ew, boolean doIMinDEE, double I0, boolean useContFlex, boolean useTupExp, boolean useEPIC, boolean checkApproxE,
        boolean outputGMECStruct, boolean eFullConfOnly, String confFileName, double stericThresh) {
        
        // TODO: refactor defaults system so we don't have to read config files to get default values
        // then we can collapse some of these ridiculously long argument lists
        // also, many settings are context-dependent
        // e.g., sometimes we want useEPIC=false and useEPIC=true in different contexts in the same Osprey run
        // see initPruning() for example
        // the all-settings-in-config-file approach doesn't makes this difficult to do
    	// in a perfect world, GMECFinder would never know what a ConfigFileParser is... we'll get there someday =)
        
        this.searchSpace = search;
        this.pruningControl = pruningControl;
        this.confSearchFactory = confSearchFactory;
        this.ecalc = ecalc;
        this.Ew = Ew;
        this.doIMinDEE = doIMinDEE;
        this.I0 = I0;
        this.useContFlex = useContFlex;
        this.useTupExp = useTupExp;
        this.useEPIC = useEPIC;
        this.checkApproxE = checkApproxE;
        this.outputGMECStruct = outputGMECStruct;
        this.EFullConfOnly = eFullConfOnly;
        this.confFileName = confFileName;
        this.stericThresh = stericThresh;
    }

    public void init(ConfigFileParser cfp) {
        init(cfp, cfp.getSearchProblem());
    }
    
    public void init(ConfigFileParser cfp, SearchProblem search) {
        
        // NOTE: this is only directly called by COMETS, which wants to use its own search problem
    	
    	// TODO: config like this would ideally be in the ConfigFileParser class
    	// so GMECFinder can live in a world where config files don't exist
        
        searchSpace = search;
        
        Ew = cfp.params.getDouble("Ew");
        doIMinDEE = cfp.params.getBool("imindee");
        if(doIMinDEE){
            I0 = cfp.params.getDouble("Ival");
        }
        
        useContFlex = cfp.params.getBool("doMinimize") || cfp.params.getBool("doPerturbations");
        //using either continuous sidechain flexibility, or backbone flexibility (which may be continuous)
        //Note discrete flexibility is just a special case of continuous flexibility
        
        useTupExp = cfp.params.getBool("UseTupExp");
        useEPIC = cfp.params.getBool("UseEPIC");
        
        checkApproxE = cfp.params.getBool("CheckApproxE");
        
        if(doIMinDEE && !useContFlex)
            throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility");
        
        outputGMECStruct = cfp.params.getBool("OUTPUTGMECSTRUCT");
        
        //for now only full-conf-only E-fcn supported is Poisson-Boltzmann
        EFullConfOnly = cfp.params.getBool("UsePoissonBoltzmann");

        usePLUG = cfp.params.getBool("UsePLUG", false);
        useDEE = cfp.params.getBool("USEDEE", true);
        if(!useDEE)
            doIMinDEE = false;//no imindee if no dee at all


        confFileName = cfp.params.getRunSpecificFileName("CONFFILENAME", ".confs.txt");
        
        // NOTE: we'll change some of these params before actually running pruning
        stericThresh = cfp.params.getDouble("StericThresh");
        
        pruningControl = new PruningControl(
            searchSpace,
            0, // pruning interval, set by initPruning()
            cfp.params.getBool("TYPEDEP"), 
            cfp.params.getDouble("BOUNDSTHRESH"),
            cfp.params.getInt("ALGOPTION"), 
            cfp.params.getBool("USEFLAGS"),
            cfp.params.getBool("USETRIPLES"),
            false,
            false, // useEPIC, set by initPruning()
            false, // useTupExp, set by initPruning()
            stericThresh
        );
        
        //initialize some kind of combinatorial search, like A*
        //FOR NOW just using A*; may also want BWM*, WCSP, or something according to settings
        confSearchFactory = ConfSearchFactory.Tools.makeFromConfig(searchSpace, cfp);
        
        // configure low-energy conformation pruning if needed
        if (cfp.params.getBool("UsePartCR")) {
            setConfPruner(new PartCRConfPruner(search, Ew));
        }
        
        // what is the "true" energy for a conformation?
	// MINIMIZED, EPIC, OR MATRIX E AS APPROPRIATE
        // TODO: these subclasses should be moved to whatever other packages care about these specific algorithms
        // TODO: let the caller set ecalc instances directly (eg in python-land)
		if (useContFlex || EFullConfOnly) {
			if ( ((useEPIC||useTupExp) && (!checkApproxE)) || searchSpace.useVoxelG ) {
				
				// use the approx minimized energy for confs
				ecalc = new GMECConfEnergyCalculator.Async.Adapter(new GMECConfEnergyCalculator() {
					@Override
					public EnergiedConf calcEnergy(ScoredConf conf) {
						double energy = searchSpace.approxMinimizedEnergy(conf.getAssignments());
						return new EnergiedConf(conf, energy);
					}
				});
				
			} else {
				
				boolean avoidCopyingMolecules = EFullConfOnly;
                                //MH: All the current conformational perturbations as of 9/12/16 should support copying
                                //to new molecules, but I'll leave this option in case new DOFs or something cause an issue
                                //Poisson-Boltzmann still an issue though...
				if (avoidCopyingMolecules) {
					System.out.println("\n\nWARNING: concurrent minimizations disabled\n");
					
					// fall back to the old SearchProblem.minimize() method
					ecalc = new GMECConfEnergyCalculator.Async.Adapter(new GMECConfEnergyCalculator() {
						@Override
						public EnergiedConf calcEnergy(ScoredConf conf) {
							double energy = searchSpace.minimizedEnergy(conf.getAssignments());
							return new EnergiedConf(conf, energy);
						}
					});
					
				} else {
				
					// for "regular" conf minimization, use the spiffy new ConfMinimizer!
					ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
					Parallelism parallelism = Parallelism.makeFromConfig(cfp);
					ecalc = MinimizingConfEnergyCalculator.make(ffparams, search, parallelism);
				}
			}
			
		} else {
			
			// for rigid calc w/ pairwise E-mtx, can just use the conf score as the energy
			ecalc = new GMECConfEnergyCalculator.Async.Adapter(new GMECConfEnergyCalculator() {
				@Override
				public EnergiedConf calcEnergy(ScoredConf conf) {
					return new EnergiedConf(conf, conf.getScore());
				}
			});
		}
    }
    
    public void setLogConfsToConsole(boolean val) {
        logConfsToConsole = val;
    }
    
    public void setConfPruner(ConfPruner val) {
        confPruner = val;
    }
    
    public List<EnergiedConf> calcGMEC(){
        return calcGMEC(I0);
    }
    
    public double calcGMECEnergy(){
        return calcGMEC().get(0).getEnergy();
    }
    
    private List<EnergiedConf> calcGMEC(double interval) {
        
        System.out.println("Calculating GMEC with interval = " + interval);
        
        boolean printEPICEnergy = checkApproxE && useEPIC && useTupExp;
        ConfPrinter confPrinter = new ConfPrinter(searchSpace, confFileName, printEPICEnergy);
        
        // 11/11/2015 JJ: This logic belongs out here. A function that does nothing if a flag is false should 
        // have its flag promoted outside of the function, unless it's used multiple times. In that case
        // the function needs to be named accordingly.
        if (useEPIC) {
            checkEPICThresh2(interval);//Make sure EPIC thresh 2 matches current interval
        }

        //precompute the energy, pruning, and maybe EPIC or tup-exp matrices
        //must be done separately for each round of iMinDEE
        precomputeMatrices(Ew + interval);
        
        // start searching for the min score conf
        System.out.println("Searching for min score conformation...");
        Stopwatch minScoreStopwatch = new Stopwatch().start();
        ConfSearch confSearch = confSearchFactory.make(searchSpace.emat, searchSpace.pruneMat);
        try {
            System.out.println("\t(among " + confSearch.getNumConformations().floatValue() + " possibilities)");
        } catch (UnsupportedOperationException ex) {
            // conf tree doesn't support it, no big deal
        }
        ScoredConf minScoreConf = confSearch.nextConf();
        if (minScoreConf == null) {
            
            // no confs in the search space, can't recover, just bail
            System.out.println("All conformations pruned. Try choosing a larger pruning interval or steric threshold.");
            return new ArrayList<>();
        }
        System.out.println("Found min score conformation in " + minScoreStopwatch.getTime(1));
        
        // evaluate the min score conf
        System.out.println("Computing energy...");
        EnergiedConf eMinScoreConf = ecalc.calcEnergy(minScoreConf);
        confPrinter.printConf(eMinScoreConf);
        System.out.println("\nMIN SCORE CONFORMATION");
        System.out.print(confPrinter.getConfReport(eMinScoreConf));
        List<EnergiedConf> econfs = new ArrayList<>();
        econfs.add(eMinScoreConf);
        
		// estimate the top of our energy window
		// this is an upper bound for now, we'll refine it as we evaluate more structures
		final EnergyRange window = new EnergyRange(eMinScoreConf.getEnergy(), Ew);
        setWindowProgress(confSearch, window);

        // enumerate all confs in order of the scores, up to the estimate of the top of the energy window
        System.out.println("Enumerating other low-scoring conformations...");
        List<ScoredConf> lowEnergyConfs = new ArrayList<>();
        Stopwatch stopwatch = new Stopwatch().start();
        lowEnergyConfs.add(minScoreConf);
        int indexToMinimizeNext = 1;
        while (true) {
            
            ScoredConf conf = confSearch.nextConf();
            if (conf == null) {
                break;
            }
            lowEnergyConfs.add(conf);
            if (conf.getScore() >= window.getMax()) {
                break;
            }
            
            // if we've been enumerating confs for a while, try a minimization to see if we get a smaller window
            if (stopwatch.getTimeS() >= 10) {
                stopwatch.stop();
                
                // save the conf and the energy for later
                EnergiedConf econf = ecalc.calcEnergy(lowEnergyConfs.get(indexToMinimizeNext++));
                handleEnergiedConf(econfs, confPrinter, window, econf);
                
                boolean changed = window.updateMin(econf.getEnergy());
                if (changed) {
                    System.out.println(String.format("Lower conformation energy updated energy window! remaining: %14.8f", window.getMax() - conf.getScore()));
                    setWindowProgress(confSearch, window);
                }
                
                stopwatch.start();
            }
        }
        System.out.println(String.format("\tFound %d more", lowEnergyConfs.size() - 1));
        
        // we're done with A*, release the tree so we can get the memory back
        confSearch = null;

        if (!lowEnergyConfs.isEmpty()) {

                // prune the confs list
                if (confPruner != null) {
                        confPruner.prune(lowEnergyConfs, ecalc);
                }

                // calculate energy for each conf
                // this will probably take a while, so track progress
                final Progress progress = new Progress(lowEnergyConfs.size());
                progress.setProgress(indexToMinimizeNext);

                // what to do when we get a conf energy?
                GMECConfEnergyCalculator.Async.Listener ecalcListener = new GMECConfEnergyCalculator.Async.Listener() {
                        @Override
                        public void onFinished(EnergiedConf econf) {

                                handleEnergiedConf(econfs, confPrinter, window, econf);
                                progress.incrementProgress();

                                // refine the estimate of the top of the energy window
                                boolean changed = window.updateMin(econf.getEnergy());
                                if (changed) {

                                        // prune conformations with the new window
                                        for (int i=lowEnergyConfs.size()-1; i>=0; i--) {
                                                if (lowEnergyConfs.get(i).getScore() > window.getMax()) {
                                                        lowEnergyConfs.remove(i);
                                                } else {
                                                        break;
                                                }
                                        }

                                        // update progress
                                        System.out.println(String.format("\nNew lowest energy: %.6f", window.getMin()));
                                        System.out.println(String.format("\tReduced to %d low-energy conformations", lowEnergyConfs.size()));
                                        progress.setTotalWork(lowEnergyConfs.size());
                                }
                        }
                };

                // calc the conf energy asynchronously
                System.out.println(String.format("\nComputing energies for %d conformations...", lowEnergyConfs.size()));
                for (; indexToMinimizeNext<lowEnergyConfs.size(); indexToMinimizeNext++) {
                        ecalc.calcEnergyAsync(lowEnergyConfs.get(indexToMinimizeNext), ecalcListener);
                }
                ecalc.getTasks().waitForFinish();
        }
		
        // sort all the confs by energy
        Collections.sort(econfs, new Comparator<EnergiedConf>() {
                @Override
                public int compare(EnergiedConf a, EnergiedConf b) {
                        return Double.compare(a.getEnergy(), b.getEnergy());
                }
        });
        
        // get the min energy conf
        EnergiedConf minEnergyConf = econfs.get(0);
        
        if(doIMinDEE){//iMinDEE...figure out if a second round is needed
            
            //compute the lowest lower-bound on a conformation
            // if we're using epic or tuple expansion, we need to compute the min bound using the energy matrix
            // otherwise, our pruning interval estimate will be wrong
            if ((useEPIC||useTupExp))//enumeration is by approximated energy...calculate lower bound separately
                lowestBound = lowestPairwiseBound(searchSpace);
            else//enumeration is by lower bound, so use that
                lowestBound = minScoreConf.getScore();
        
            // could the minGMEC have been pruned due to a pruning interval that's too small?
            if (minEnergyConf.getEnergy() > lowestBound + interval) {

                // yeah, it could have been. we can't prove minEnergyConf is the minGMEC
                // we have to pick a new interval and try again
                System.out.println("Pruning interval is too small. minGMEC could have been pruned.");
                System.out.println("Will estimate new interval based on conformations evaluated so far and restart");


                double nextInterval = minEnergyConf.getEnergy() - lowestBound;

                // pad the new interval a bit to avoid numerical instability
                double intervalPad = 0.001;
                if(searchSpace.useVoxelG)//GMEC energy has statistical error (~0.1 kcal/mol maybe)
                    intervalPad = 0.2;
                nextInterval += intervalPad;

                return calcGMEC(nextInterval);
            }
        }
        
        // we didn't prune it! minEnergyConf is the minGMEC!! =)
        EnergiedConf minGMEC = minEnergyConf;
        System.out.println(String.format("\nFound minGMEC!"));
        System.out.print(confPrinter.getConfReport(minGMEC));
        if (outputGMECStruct) {
            searchSpace.outputMinimizedStruct(minGMEC.getAssignments(), searchSpace.name + ".GMEC.pdb");
        }
        
        // prune all confs outside the energy window and return them
        Iterator<EnergiedConf> iter = econfs.iterator();
        while (iter.hasNext()) {
            EnergiedConf econf = iter.next();
            if (econf.getEnergy() > minGMEC.getEnergy() + Ew) {
                iter.remove();
            }
        }
        
        if (Ew > 0) {
            System.out.println(String.format("Also found %d more conformations in energy window", econfs.size() - 1));
        }
        
        if (ecalc instanceof MinimizingConfEnergyCalculator) {
        	((MinimizingConfEnergyCalculator)ecalc).clean();//Clean up once both rounds of iMinDEE (if applicable) are done
        }
        return econfs;
    }
    
    private void handleEnergiedConf(List<EnergiedConf> econfs, ConfPrinter confPrinter, EnergyRange window, EnergiedConf econf) {
        
        econfs.add(econf);
        
        // immediately output the conf, in case the run aborts and we want to resume later
        confPrinter.printConf(econf);

        // log the conf to console if desired
        if (logConfsToConsole) {
            System.out.println("\nENUMERATING CONFORMATION");
            System.out.print(confPrinter.getConfReport(econf, window));
        }
    }

    private void setWindowProgress(ConfSearch confSearch, EnergyRange window) {
        
        // HACKHACK: set progress goal
        if (confSearch instanceof ConfAStarTree) {
            ConfAStarTree tree = (ConfAStarTree)confSearch;
            if (tree.getProgress() != null) {
                tree.getProgress().setGoalScore(window.getMax());
            }
        }
    }

    public List<ScoredConf> calcSequences() {
        return calcSequences(I0);
    }
    
    public List<ScoredConf> calcSequences(double interval) {
        
        System.out.println("Finding sequences with interval " + interval + " and energy window " + Ew + " ...");
        
        boolean printEPICEnergy = checkApproxE && useEPIC && useTupExp;
        ConfPrinter confPrinter = new ConfPrinter(searchSpace, confFileName, printEPICEnergy);
        
        if (useEPIC) {
            checkEPICThresh2(interval);//Make sure EPIC thresh 2 matches current interval
        }
        
        // precompute the energy, pruning, and maybe EPIC or tup-exp matrices
        precomputeMatrices(Ew + interval);
        
        // start searching for the min score conf
        System.out.println("Searching for min score conformation...");
        Stopwatch minScoreStopwatch = new Stopwatch().start();
        ConfSearch confSearch = confSearchFactory.make(searchSpace.emat, searchSpace.pruneMat);
        try {
            System.out.println("\t(among " + confSearch.getNumConformations().floatValue() + " possibilities)");
        } catch (UnsupportedOperationException ex) {
            // conf tree doesn't support it, no big deal
        }
        ScoredConf minScoreConf = confSearch.nextConf();
        if (minScoreConf == null) {
            
            // no confs in the search space, can't recover, just bail
            System.out.println("All conformations pruned. Try choosing a larger pruning interval.");
            return new ArrayList<>();
        }
        System.out.println("Found min score conformation in " + minScoreStopwatch.getTime(1));
        
        // evaluate the min score conf
        System.out.println("Computing energy...");
        EnergiedConf eMinScoreConf = ecalc.calcEnergy(minScoreConf);
        
        // start the sequence list with the min score conf
        System.out.println("\nMIN SCORE CONFORMATION");
        System.out.print(confPrinter.getConfReport(eMinScoreConf));
        List<ScoredConf> sequenceConfs = new ArrayList<>();
        Set<String> sequenceKeys = new HashSet<>();
        sequenceConfs.add(minScoreConf);
        sequenceKeys.add(makeSequenceKey(minScoreConf));
        
        // estimate the top of our energy window
        // this is an upper bound for now, we'll refine it as we evaluate more structures
        final EnergyRange window = new EnergyRange(eMinScoreConf.getEnergy(), Ew);
        
        // enumerate all confs in order of the scores, up to the estimate of the top of the energy window
        System.out.println("Enumerating other low-scoring conformations...");
        while (true) {
            
            ScoredConf conf = confSearch.nextConf();
            if (conf == null) {
                break;
            }
            
            // is this a new sequence?
            boolean isNewSequence = sequenceKeys.add(makeSequenceKey(conf));
            if (isNewSequence) {
                
                // record the new sequence
				if (logConfsToConsole) {
					System.out.println("\nUNIQUE SEQUENCE " + sequenceConfs.size());
					System.out.print(confPrinter.getConfReport(conf, window));
				}
                sequenceConfs.add(conf);
            }
            
            if (conf.getScore() >= window.getMax()) {
                break;
            }
        }
        
        return sequenceConfs;
    }
    
    private String makeSequenceKey(ScoredConf conf) {
        StringBuilder buf = new StringBuilder();
        
        for (int pos=0; pos<searchSpace.confSpace.numPos; pos++) {
            
            String aaType = searchSpace.confSpace.posFlex.get(pos).RCs.get(conf.getAssignments()[pos]).AAType;
            
            if (buf.length() > 0) {
                buf.append(",");
            }
            buf.append(aaType);
        }
        
        return buf.toString();
    }

    //for use after an iMinDEE run, if we want to see what the lowest bound was
    public double getLowestBound(){
        return lowestBound;
    }
    
    
    public void precomputeMatrices(double pruningInterval){
        //Precalculate TupleMatrices needed for GMEC computation.  Some of these may already be computed.  
        //All of these matrices except the basic pairwise energy matrix are pruning-dependent:
        //we can prune conformations whose energies are within pruningInterval
        //of the lowest pairwise lower bound
        
        if(EFullConfOnly){//Perform a tuple expansion that does not assume a pairwise
            //energy function, and thus must omit some of the usual pruning steps.
            fullConfOnlyTupExp();
            return;
        }
        else if(searchSpace.useVoxelG){//Can't do DEE here either, but want EPIC
            voxelGTupExp();
            return;
        }
    
        //First calculate the pairwise energy matrix, if not already present
        if (searchSpace.emat == null) {
            searchSpace.loadEnergyMatrix();
        }
        
        
        
        //DEBUG!!! Timing/profiling minimization for 40.cont...about 1/3 SAPE, small EPIC terms probably unimportant...
        /*searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace,0);
        searchSpace.loadEPICMatrix();
        long startTime = System.currentTimeMillis();
        int conf[] = new int[] { 0, 24, 1, 6, 4, 0, 9, 11, 4, 4, 4, 9, 3, 2, 22, 6, 
            4, 1, 4, 1, 4, 6, 13, 0, 3, 4, 1, 3, 1, 3, 2, 6, 21, 8, 5, 0, 26, 3, 1, 19 };

        //int conf[] = new int[] {0, 26, 1, 6, 2, 0, 7, 8, 6, 4, 4, 14, 4, 3, 26, 6, 3, 1, 2, 1, 
        //    5, 9, 13, 1, 2, 4, 1, 3, 2, 7, 2, 6, 17, 21, 4, 0, 15, 6, 1, 32 };
        
        
        System.out.println("EPIC MIN: "+searchSpace.EPICMinimizedEnergy(conf));
        long minDoneTime = System.currentTimeMillis();
        System.out.println("MIN TIME: "+(minDoneTime-startTime));
        System.exit(0);*/
        

        //Doing competitor pruning now
        //will limit us to a smaller, but effective, set of competitors in all future DEE
        if(searchSpace.competitorPruneMat == null && useDEE){
            System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
            initPruning(0, false, false);

            if(usePLUG){
                searchSpace.loadPLUGMatrix();
                searchSpace.plugMat.doMultiTermPruning(searchSpace.pruneMat, PLUGPruneTriples);
            }

            pruningControl.setOnlyGoldstein(true);
            pruningControl.prune();
            searchSpace.competitorPruneMat = searchSpace.pruneMat;
            searchSpace.pruneMat = null;
            System.out.println("COMPETITOR PRUNING DONE");
        }
        
        
        //Next, do DEE, which will fill in the pruning matrix
        initPruning(pruningInterval, false, false);
        if(usePLUG){//DEBUG!!  should reuse from competitor prolly
            searchSpace.loadPLUGMatrix();
            searchSpace.plugMat.doMultiTermPruning(searchSpace.pruneMat, PLUGPruneTriples);
        }

        if(useDEE)
            pruningControl.prune();//pass in DEE options, and run the specified types of DEE
        
        
        //precomputing EPIC or tuple-expander matrices is much faster
        //if only done for unpruned RCs.  Less RCs to handle, and the fits are far simpler.  
        if(useEPIC){
            searchSpace.loadEPICMatrix();
            
            //we can prune more using the EPIC matrix
            if(searchSpace.epicSettings.useEPICPruning && useDEE){
                System.out.println("Beginning post-EPIC pruning.");
                initPruning(pruningInterval, true, false);
                pruningControl.prune();
                System.out.println("Finished post-EPIC pruning.");
            }
        }
        if(useTupExp){//preferably do this one EPIC loaded (much faster if can fit to EPIC)
            searchSpace.loadTupExpEMatrix();

            if(useDEE) {
                //we can prune even more with tup-exp!
                //we can prune more using the EPIC matrix
                //no iMinDEE interval needed here
                System.out.println("Beginning post-tup-exp pruning.");
                initPruning(Ew, false, true);
                pruningControl.prune();
                System.out.println("Finished post-tup-exp pruning.");
            }
        }
    }
    
 
    private void initPruning(double pruningInterval, boolean useEPIC, boolean useTupExp) {
        
        // init the pruning matrix if needed
        if(searchSpace.pruneMat == null || searchSpace.pruneMat.getPruningInterval() < pruningInterval) {
            searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, pruningInterval);
        }
        
        // configure the pruner
        pruningControl.setOnlyGoldstein(false);
        pruningControl.setPruningInterval(pruningInterval);
        pruningControl.setUseEPIC(useEPIC);
        pruningControl.setUseTupExp(useTupExp);
    }
 
    private void fullConfOnlyTupExp(){
        //precompute the tuple expansion
        if(!useTupExp)
            throw new RuntimeException("ERROR: Need tuple expansion to handle full-conf-only E-function");
        if(useEPIC)//later consider using differencing scheme to do EPIC for these
            throw new RuntimeException("ERROR: EPIC for full-conf-only E-function not yet supported");
        if(doIMinDEE)//don't have concept of pairwise lower bound, so not doing iMinDEE 
            //(can just do rigid pruning on tup-exp matrix, even if using cont flex)
            throw new RuntimeException("ERROR: iMinDEE + full-conf-only E-function not supported");
        
        
        //Let's compute a matrix from the pairwise terms (no P-B), to use in selecting triples
        searchSpace.loadEnergyMatrix();
        
        //initialize pruning matrix.  Nothing pruned yet because don't have pairwise energies
        searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace,Ew);//not iMinDEE
        
        //We can only do steric pruning
        //May want to set a lower thresh than the default (30 perhaps)
        Pruner pruner = new Pruner(searchSpace, false, 0, 0, false, false);
        pruner.pruneSteric(stericThresh);

        if(usePLUG){//still valid here, and provides better pruning
            searchSpace.loadPLUGMatrix();
            searchSpace.plugMat.doMultiTermPruning(searchSpace.pruneMat, PLUGPruneTriples);
        }
                
        searchSpace.loadTupExpEMatrix();
    }


    private void voxelGTupExp(){
        //precompute the tuple expansion
        if(!useTupExp)
            throw new RuntimeException("ERROR: Need tuple expansion to handle voxelG");
        if(!useEPIC)//later consider using differencing scheme to do EPIC for these
            throw new RuntimeException("ERROR: Need EPIC to handle voxelG");
        if(doIMinDEE)//don't have concept of pairwise lower bound, so not doing iMinDEE
            //(can just do rigid pruning on tup-exp matrix, even if using cont flex)
            throw new RuntimeException("ERROR: iMinDEE + voxelG not supported");


        //Let's compute a matrix from the pairwise terms (no entropy), to use in selecting triples
        searchSpace.loadEnergyMatrix();

        //initialize pruning matrix.  Nothing pruned yet because don't have pairwise energies
        searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace,Ew);//not iMinDEE

        Pruner pruner = new Pruner(searchSpace, false, 0, 0, false, false);
        pruner.pruneSteric(stericThresh);

        if(usePLUG){//still valid here, and provides better pruning
            searchSpace.loadPLUGMatrix();
            searchSpace.plugMat.doMultiTermPruning(searchSpace.pruneMat, PLUGPruneTriples);
        }

        searchSpace.loadEPICMatrix();//sets up both EPIC matrix and voxel G calculator
        searchSpace.loadTupExpEMatrix();
    }

    
    public double lowestPairwiseBound(SearchProblem prob){
        //In an EPIC calculation, our enumeration will probably include much less conformations,
        //but for iMinDEE purposes we still need to know what our lowest bound would have been
        //if we enumerated w/o EPIC (i.e., what is the minimum energy calculated using the lower-bound energy matrix)
        
        System.out.println();
        System.out.println("Calculating no-EPIC lowest energy bound");
        
        SearchProblem probNoEPIC = new SearchProblem(prob);//shallow copy
        probNoEPIC.useEPIC = false;
        probNoEPIC.useTupExpForSearch = false;
        
        ConfSearch searchNoEPIC = ConfTree.makeFull(probNoEPIC);
        ScoredConf LBConf = searchNoEPIC.nextConf();//lowest conf for this search
        
        double LB = Double.POSITIVE_INFINITY;
        if(LBConf != null)
            LB = probNoEPIC.lowerBound(LBConf.getAssignments());
        
        System.out.println("No-EPIC lowest energy bound: "+LB);
        System.out.println();
        
        return LB;
    }

    
    private void checkEPICThresh2(double curInterval){
        if(curInterval+Ew>searchSpace.epicSettings.EPICThresh2){//need to raise EPICThresh2 
            //to the point that we can guarantee no continuous component of the GMEC
            //or desired ensemble will need to reach it
            System.out.println("Raising EPICThresh2 to "+(curInterval+Ew)+" based on "
                    + "iMinDEE interval and energy window");
            searchSpace.epicSettings.EPICThresh2 = curInterval+Ew;
        }
    }
   
}
