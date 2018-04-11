package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;

/**
 *
 * @author lowegard (adapted from GMECFinder.java)
 */

public class EWAKStarLowEnergyFinder {

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

    private double I0=0;//initial value of iMinDEE pruning interval
    private boolean doIMinDEE;//do iMinDEE

    private boolean useContFlex;
    //boolean enumByPairwiseLowerBound;//are we using a search method that
    //enumerates by pairwise lower bound (minDEE-style)
    //or by (possibly approximated) true energy (rigid, EPIC, or tuple expander)?

    private double stericThresh;


    public EWAKStarLowEnergyFinder() {

        // Arguably, the stuff below is the initialization, which should be its own function. The constructor may eventually
        // do other more interesting things and reading in member variables isn't that interesting. -JJ
        // Ew = cfp.prams.getDouble("Ew");
        // etc...

        // Jeff: I think is a good place to set default values
        ecalc = null;

        // TODO: in the future, we want all config options to act like the ones set in this constructor
        // The goal is: construct a GMECFinder instance without huge arguments lists and it will use defaults automatically
        // then, when you want to deviate from defaults, you call e.g. :  gmecFinder.setThisOption(toThisVal)
        // this will make the eventual python scripting interface much much nicer!
        // options that can't fundamentally be defaulted are actually inputs, and they should generally get passed
        // in the method args in the method that actually needs them
    }

    public void init(ConfigFileParser cfp) {

        init(cfp, cfp.getSearchProblem());
    }

    public void init(ConfigFileParser cfp, SearchProblem search) {

        // NOTE: this is only directly called by COMETS, which wants to use its own search problem

        // TODO: config like this would ideally be in the ConfigFileParser class
        // so GMECFinder can live in a world where config files don't exist

        searchSpace = search;
        doIMinDEE = cfp.params.getBool("imindee");
        if(doIMinDEE){
            I0 = cfp.params.getDouble("Ival");
        }

        useContFlex = cfp.params.getBool("doMinimize") || cfp.params.getBool("doPerturbations");
        //using either continuous sidechain flexibility, or backbone flexibility (which may be continuous)
        //Note discrete flexibility is just a special case of continuous flexibility

        if(doIMinDEE && !useContFlex)
            throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility");

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

        // for "regular" conf minimization, use the spiffy new ConfMinimizer!
        ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
        Parallelism parallelism = Parallelism.makeFromConfig(cfp);
        ecalc = MinimizingConfEnergyCalculator.make(ffparams, search, parallelism);

    }

    public EnergiedConf calcSequences() {
        return calcSequences(I0);
    }

    private EnergiedConf calcSequences(double interval) {

        System.out.println("Finding a low energy conformation for the current sequence.");

        // precompute the energy, pruning, and maybe EPIC or tup-exp matrices
        precomputeMatrices(interval);

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
            return new EnergiedConf(null, 0);
        }

        System.out.println("Found min score conformation in " + minScoreStopwatch.getTime(1));

        // evaluate the min score conf
        System.out.println("Computing energy...");
        EnergiedConf eMinScoreConf = ecalc.calcEnergy(minScoreConf);

        if (eMinScoreConf.getEnergy() >= 0){
            System.out.println("The first conformation has an energy greater than zero. Continuing search... ");
        }

        //check to see if the sequence you've found has a negative energy
        //If you haven't, keep enumerating conformations and keep the first conformation with negative energy
        while (eMinScoreConf.getEnergy() >= 0) {

            ScoredConf conf = confSearch.nextConf();
            if (conf == null) {
                break;
            }

            eMinScoreConf = ecalc.calcEnergy(conf);

        }

        // start the sequence list with the min score conf
        System.out.println("\nFOUND A LOW-ENERGY CONFORMATION.");
        System.out.print(printConf(eMinScoreConf));

        return eMinScoreConf;

    }

    private String printConf(EnergiedConf conf){

        StringBuilder buf = new StringBuilder();
        int LabelSize = 30;
        String LabelFormat = "\t%-" + LabelSize + "s";

        buf.append("RESTYPES: ");
        for (int pos=0; pos<searchSpace.confSpace.numPos; pos++) {
            String resType = searchSpace.confSpace.posFlex.get(pos).RCs.get(conf.getAssignments()[pos]).AAType;
            buf.append(String.format(" %3s", resType));
        }

        buf.append(String.format(LabelFormat + " %.6f", "Energy:", conf.getEnergy()));
        buf.append(String.format(LabelFormat + " %.6f (gap: %.6f", "Score", conf.getScore(), Math.abs(conf.getScore() - conf.getEnergy())));

        buf.append(")\n");

        return buf.toString();
    }


    private void precomputeMatrices(double pruningInterval){
        //Precalculate TupleMatrices needed for GMEC computation.  Some of these may already be computed.
        //All of these matrices except the basic pairwise energy matrix are pruning-dependent:
        //we can prune conformations whose energies are within pruningInterval
        //of the lowest pairwise lower bound


        //First calculate the pairwise energy matrix, if not already present
        if (searchSpace.emat == null) {
            searchSpace.loadEnergyMatrix();
        }

        //Doing competitor pruning now
        //will limit us to a smaller, but effective, set of competitors in all future DEE
        if(searchSpace.competitorPruneMat == null){
            System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
            initPruning(0);
            pruningControl.setOnlyGoldstein(true);
            pruningControl.prune();
            searchSpace.competitorPruneMat = searchSpace.pruneMat;
            searchSpace.pruneMat = null;
            System.out.println("COMPETITOR PRUNING DONE");
        }

        //Next, do DEE, which will fill in the pruning matrix
        initPruning(pruningInterval);
        pruningControl.prune();//pass in DEE options, and run the specified types of DEE

    }

    private void initPruning(double pruningInterval) {

        // init the pruning matrix if needed
        if(searchSpace.pruneMat == null || searchSpace.pruneMat.getPruningInterval() < pruningInterval) {
            searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, pruningInterval);
        }

        // configure the pruner
        pruningControl.setOnlyGoldstein(false);
        pruningControl.setPruningInterval(pruningInterval);
        pruningControl.setUseEPIC(false);
        pruningControl.setUseTupExp(false);
    }

}