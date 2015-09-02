/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

/**
 *
 * @author mhall44
 */
class COMETSDoer {
    
    /*
    ArrayList<SearchProblem> stateSP;//search problems for each state
    
    public COMETSDoer (ConfigFileParser cfgP){
        //fill in all the settings
        
        cfp = cfgP;
        
        Ew = cfgP.params.getDouble("Ew",0);
        doIMinDEE = cfgP.params.getBool("imindee",false);
        if(doIMinDEE){
            I0 = cfgP.params.getDouble("Ival",5);
        }
        
        useContFlex = cfgP.params.getBool("doMinimize",false);
        useTupExp = cfgP.params.getBool("UseTupExp",false);
        useEPIC = cfgP.params.getBool("UseEPIC",false);
        
        checkApproxE = cfgP.params.getBool("CheckApproxE",true);
        
        if(doIMinDEE && !useContFlex)
            throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility");
        
        outputGMECStruct = cfgP.params.getBool("OUTPUTGMECSTRUCT", false);
        
        //for now only full-conf-only E-fcn supported is Poisson-Boltzmann
        EFullConfOnly = cfgP.params.getBool("UsePoissonBoltzmann",false);
    }
    
    public void calcBestSequences(){
        see GMECFinder;
        //OK if cont flex with LUTE i-vals would be good
        //let's start discrete though
        //with continuous i guess each E_a(s) technically has its own ival...
        //we want the max one...
        /*
        Yeah...the IVal we use for each state must be the upper bound on sequence Ivals for that state
        So, it is valid to use any upper bound on max opt-sequence GMEC for that state
        - lowestBound for that (over all sequences)
        this is fairly easy if we have stability constraints
        Another option is to just use local-clash-avoidance pruning
        Also could consider making a mapping from bound to unbound E-terms
        and bound that to figure out how much worse unbound can be compared to bound...
        Also we can deem seqs to be clashing above particular ival...
        finally we could consider just doing minDEE and then LUTE.  Valid too
        and doesn't require cont flex in EPIC itself
        */
    //}
    
}
