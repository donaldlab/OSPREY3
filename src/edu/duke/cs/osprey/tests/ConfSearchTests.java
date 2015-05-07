/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchSpace;
import edu.duke.cs.osprey.pruning.PruningControl;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class ConfSearchTests {
    //Making sure each type of ConfSearch minimizes properly, 
    //and that DEE pruning does not change the overall minimum
    //this is testing discrete conformational search
    
    public static void testExhaustive(){
        //For a small, simple search space, ensure that the ConfSearch object
        //enumerates all conformations in ascending order
        
        SearchSpace searchSpace = makeTestSearchSpace(5,false,false);//3^5 confs.  No clashes (would make some weird confs)
        int totNumResults = 3*3*3*3*3;
        
        ArrayList<ConfSearch> searches = confSearchesToTest(searchSpace);
        
        for(ConfSearch search : searches){
            int resultCount=0;
            
            double curE = Double.NEGATIVE_INFINITY;
            
            while(true){
                int[] conf = search.nextConf();
                if(conf==null)
                    break;
                
                double newConfE = searchSpace.emat.confE(conf);
                
                assert newConfE > curE - 1e-3;//enumerating in ascending order (within numerical tolerance)
                curE = Math.max(curE,newConfE);
            }
            
            assert resultCount==totNumResults;//check right number of confs enumerated
        }
        
        System.out.println("EXHAUSTIVE CONFORMATIONAL SEARCH TEST PASSED");
    }
    
    
    public static void testDEE(){
        //ensure that DEE does not prune the best conf, by comparing A* results with and without DEE
        //might run this a few time with random energies, to be more sure
        
        SearchSpace searchSpace = makeTestSearchSpace(6/*10*/,true,false);
        //bigger search space, and leaving clashes will create more realistic test conditions
        
        ConfSearch aStar = new ConfTree(searchSpace);//Regular A* is cool for this purpose
        
        int topConf[] = aStar.nextConf();
        
        int algOption = 1;
        double boundsThresh = Double.POSITIVE_INFINITY;//could also try energy of topConf...
        
        //now prune and rerun A*
        PruningControl pruning = new PruningControl(searchSpace, 0, false, boundsThresh,
                algOption, true, false, false);
                
        pruning.prune();
        
        
        int topConfWithPruning[] = aStar.nextConf();
        
        for(int pos=0; pos<searchSpace.confSpace.numPos; pos++)
            assert topConf[pos] == topConfWithPruning[pos];
        
        System.out.println("DEE TEST PASSED");
    }
    
    private static ArrayList<ConfSearch> confSearchesToTest(SearchSpace searchSpace){
        //for a given search space, enumerate the ConfSearch objects we want to test for it
        ArrayList<ConfSearch> ans = new ArrayList<>();
        ans.add(new ConfTree(searchSpace));
        //BWM*, WCSP, etc. here?  Only provable methods to be considered
        return ans;
    }
    
    
    
    
    private static SearchSpace makeTestSearchSpace(int numPos, boolean doMut, boolean randomizeEnergies){
        //generate a search space for test purposes with the given number of positions
        //If doMut, allow a bunch of options per position, else allow only val (3 RCs) at each position
        //We'll randomize the energies to get a fresh test every time, but 
        //if allowClashes the really big energies will be left alone (to simulate
        //real energies that have clashes)
        
        //options for mutations at each position
        ArrayList<String> AAatPos = new ArrayList<>();
        AAatPos.add("Val");//3 rotamers
        if(doMut){//34 rotamers in all
            AAatPos.add("Phe");
            AAatPos.add("Lys");
        }
        
        ArrayList<String> flexRes = new ArrayList<>();
        ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
        
        for(int pos=0; pos<numPos; pos++){
            //we'll start with residue 20 and go up.  We assume we won't hit the end of the molecule
            flexRes.add(Integer.valueOf(20+pos).toString());
            allowedAAs.add(AAatPos);
        }
        
        SearchSpace ans = new SearchSpace( "CONFSEARCHTEST", "1CC8.ss.pdb", 
                flexRes, allowedAAs,
                false, false );//don't add WT, and no minimization
        
        if(randomizeEnergies){
            //we don't need real energies, just make some up (in fact the randomization will be good)
            for(ArrayList<Double> resE : ans.emat.oneBody)
                fillInRandomly(resE);
            for(ArrayList<ArrayList<ArrayList<Double>>> resE : ans.emat.pairwise){
                for(ArrayList<ArrayList<Double>> pairE : resE){
                    for(ArrayList<Double> rotE : pairE)
                        fillInRandomly(rotE);
                }
            }
        }
        else {
            ans.loadEnergyMatrix();
        }
        
        return ans;
    }
    
    
    private static void fillInRandomly(ArrayList<Double> energies){
        //randomly fill in a list of "energies".  Keep them within a reasonable range
        for(int index=0; index<energies.size(); index++){
            energies.set( index, 20*(Math.random()-0.5) );
        }
    }
    
    
    
}
