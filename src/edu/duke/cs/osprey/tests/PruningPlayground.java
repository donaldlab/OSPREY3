/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.TupleEnumerator;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * the main game here is comparing PLUG and DEE
 *
 * @author mhall44
 */
public class PruningPlayground {

    public static void main(String args[]){
        Switches sw = new Switches(args);
        if(sw.usingSwitches())
            args = Arrays.copyOfRange(args, 0, args.length-1);
        PruningPlayground fun = new PruningPlayground(args);
        if(sw.isSwitchedOn(-1,false)){
            System.out.println("Single-term PLUG only: ");
            fun.checkPruningEfficacy(true, false, false, 0);
        }
        if(sw.isSwitchedOn(-1,true)){
            System.out.println("PLUG only: ");
            fun.checkPruningEfficacy(true, true, false, 0);
        }
        for(double ival : new double[] {0,5,10,15,20,1000000}){//top one is to just do steric
            System.out.println("Ival: "+ival);
            if(sw.isSwitchedOn(ival,false)){
                System.out.println("DEE only: ");
                fun.checkPruningEfficacy(false, false, true, ival);
            }
            if(sw.isSwitchedOn(ival,true)){
                System.out.println("DEE+PLUG: ");
                fun.checkPruningEfficacy(true, true, true, ival);
            }
        }
    }


    String args[];
    PolytopeMatrix plugMatCache = null;

    public PruningPlayground(String[] args){
        this.args = args;
    }


    private void checkPruningEfficacy(boolean usePLUG, boolean fullPLUGPruning, boolean useDEE, double pruningInterval){

        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(args);//args 1, 3+ are configuration files
        cfp.loadData();


        SearchProblem sp = cfp.getSearchProblem();


        sp.loadEnergyMatrix();
        sp.pruneMat = new PruningMatrix(sp.confSpace, pruningInterval);

        //if doing both PLUG and DEE, do PLUG first because it can help DEE
        if(usePLUG){
            if(fullPLUGPruning){
                if(plugMatCache==null)
                    plugMatCache = new PolytopeMatrix(sp, true);
                plugMatCache.doMultiTermPruning(sp.pruneMat, true);//includes all single-term pruning
            }
            else
                new PolytopeMatrix(sp, true);//don't cache or do multi-term pruning
        }
        if(useDEE){
            boolean useTriples = pruningInterval < 1000;//DEBUG!!!  basically do this unless only steric
            if(useDEE)//do this now bc if allow PLUG-pruned competitors they may prune everything (between the conflicting PLUG and regular criteria)
                prepareCompetitorPruneMat(sp, pruningInterval, useTriples);
            PruningControl pruningControl = new PruningControl(sp, pruningInterval, false, 100, 3, true, useTriples, false, false, false, 100);
            pruningControl.prune();
        }

        TupleEnumerator tupEnum = new TupleEnumerator(sp.pruneMat, null, sp.confSpace.numPos);
        System.out.println(tupEnum.enumerateUnprunedTuples(1).size()+" unpruned singles, out of "+totalNumSingles(sp));
        System.out.println(tupEnum.enumerateUnprunedTuples(2).size()+" unpruned pairs");
        /*try {
            System.out.println("Is expected GMEC pruned?: "+sp.pruneMat.isPruned(new RCTuple(new int[] {5,7,7,5,1,7,4})));
        }
        catch(Exception e) {
            System.out.println("lol not basic 1CC8");
        }*/
    }


    private static int totalNumSingles(SearchProblem sp){
        int ans = 0;
        for(int nr : sp.confSpace.getNumRCsAtPos())
            ans += nr;
        return ans;
    }


    private void prepareCompetitorPruneMat(SearchProblem sp, double pruningInterval, boolean useTriples){
        System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
        PruningMatrix tmpPruneMat = sp.pruneMat;
        sp.pruneMat = new PruningMatrix(sp.confSpace, 0);
        //want to remove everything that PLUG pruned though  DEBUG!!!
        TupleEnumerator tupEnum = new TupleEnumerator(sp.pruneMat,  null, sp.confSpace.numPos);
        for(int numPos=1; numPos<=3; numPos++){
            for(RCTuple tup : tupEnum.enumerateUnprunedTuples(numPos)){
                if(tmpPruneMat.isPruned(tup))
                    sp.pruneMat.setTupleValue(tup,true);
            }
        }


        PruningControl pruningControl = new PruningControl(sp, 0., false, 100, 3, true, useTriples, false, false, false, 100);
        pruningControl.setOnlyGoldstein(true);
        pruningControl.prune();
        sp.competitorPruneMat = sp.pruneMat;
        sp.pruneMat = tmpPruneMat;
        System.out.println("COMPETITOR PRUNING DONE");
    }



    //switches so we can run a subset of the runs (due to time limit)
    private static class Switches {
        //if the last argument of args starts with SWITCHES, interpret it as a list
        //of what conditions (ival, plug or not) to use
        //format: SWITCHES,-1,true,5,false etc (-1 meaning PLUG only)
        ArrayList<RunPair> switchedOn = null;

        private Switches(String args[]){//we expect at least one arg...there must be config files
            String arg = args[args.length-1];
            if(arg.startsWith("SWITCHES")){
                switchedOn = new ArrayList<>();
                String tok[] = arg.split(",");
                for(int a=0; a<tok.length/2; a++){
                    switchedOn.add(new RunPair(Double.parseDouble(tok[2*a+1]),Boolean.parseBoolean(tok[2*a+2])));
                }
            }
        }

        boolean usingSwitches(){
            return switchedOn!=null;
        }

        boolean isSwitchedOn(double ival, boolean second){
            if(!usingSwitches())//everything on by default
                return true;
            for(RunPair rp : switchedOn){
                if(rp.matches(ival, second))
                    return true;
            }
            return false;
        }

        //Let's have a class within a class within a class
        private static class RunPair {
            double ival;
            boolean second;//indicates the second run in the pair (e.g. DEE+PLUG)
            private RunPair(double ival, boolean second){
                this.ival = ival;
                this.second = second;
            }
            private boolean matches(double ival, boolean second){
                return Math.abs(ival-this.ival)<1e-10 && second==this.second;
            }
        }
    }



}
