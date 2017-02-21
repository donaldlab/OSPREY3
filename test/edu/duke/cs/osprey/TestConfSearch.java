/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.InfiniteIterator;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import java.util.ArrayList;
import java.util.Random;

/**
 *
 * @author mhall44
 */
public class TestConfSearch extends TestBase {
    //Checking accuracy of conformational search under various conditions,
    //and that DEE pruning does not change the overall minimum
    //this is testing discrete conformational search
    
    @BeforeClass
    public static void before() {
            initDefaultEnvironment();
    }
    
    
    @Test
    public void runTests(){
        testDEE(true);
        testDEE(false);
        testExhaustive(false, false);
        testExhaustive(false, true);
        testExhaustive(true, false);
    }
    
    public void testExhaustive(boolean useTriples, boolean useEllipses){
        //For a small, simple search space, ensure that the ConfSearch object
        //enumerates all conformations in ascending order
        
        SearchProblem searchSpace = makeTestSearchSpace(5,false,true,useTriples,useEllipses);//3^5 confs.  No clashes (would make some weird confs)
        int totNumResults = 3*3*3*3*3;
        
        searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace,-1);//no pruning
        ArrayList<ConfSearch> searches = confSearchesToTest(searchSpace);
        
        for(ConfSearch search : searches){
            int resultCount=0;
            
            double curE = Double.NEGATIVE_INFINITY;
            
            while(true){
                ConfSearch.ScoredConf conf = search.nextConf();
                if(conf==null)
                    break;
                
                assert conf.getScore() > curE - 1e-3;//enumerating in ascending order (within numerical tolerance)
                curE = Math.max(curE,conf.getScore());
                resultCount++;
            }
            
            assert (resultCount==totNumResults);//check right number of confs enumerated
        }
        
        System.out.println("EXHAUSTIVE CONFORMATIONAL SEARCH TEST PASSED.  useTriples: "+useTriples+" useEllipses: "+useEllipses);
    }
    
    
    public static void testDEE(boolean useEllipses){
        //ensure that DEE does not prune the best conf, by comparing A* results with and without DEE
        //might run this a few time with random energies, to be more sure
        
        SearchProblem searchSpace = makeTestSearchSpace(6/*10*/,true,false,false,useEllipses);

        //bigger search space, and leaving clashes will create more realistic test conditions
        
        searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace,-1);//no pruning
        ConfSearch aStar = ConfTree.makeFull(searchSpace);//Regular A* is cool for this purpose
        
        ConfSearch.ScoredConf topConf = aStar.nextConf();
        
        System.out.println("Conf E: "+topConf.getScore());
        
        int algOption = 3;
        double boundsThresh = Double.POSITIVE_INFINITY;//could also try energy of topConf...
        
        //now prune and rerun A*
        PruningControl pruning = new PruningControl(searchSpace, 0, false, boundsThresh,
                algOption, true, true, false, false, false, Double.POSITIVE_INFINITY);
                
        pruning.prune();
        
        aStar = ConfTree.makeFull(searchSpace);
        
        ConfSearch.ScoredConf topConfWithPruning = aStar.nextConf();
        System.out.println("Conf E: "+topConfWithPruning.getScore());

        
        for(int pos=0; pos<searchSpace.confSpace.numPos; pos++){
            boolean match = ( topConf.getAssignments()[pos] == topConfWithPruning.getAssignments()[pos] );
            assert match;
        }
        
        System.out.println("DEE TEST PASSED");
    }
    
    private static ArrayList<ConfSearch> confSearchesToTest(SearchProblem searchSpace){
        //for a given search space, enumerate the ConfSearch objects we want to test for it
        ArrayList<ConfSearch> ans = new ArrayList<>();
        ans.add(ConfTree.makeFull(searchSpace));
        //BWM*, WCSP, etc. here?  Only provable methods to be considered
        return ans;
    }
    
    
    
    
    private static SearchProblem makeTestSearchSpace(int numPos, boolean doMut, boolean randomizeEnergies, boolean includeTriples,
                boolean useEllipses){

        //generate a search space for test purposes with the given number of positions
        //If doMut, allow a bunch of options per position, else allow only val (3 RCs) at each position
        //We'll randomize the energies to get a fresh test every time, but 
        //if allowClashes the really big energies will be left alone (to simulate
        //real energies that have clashes)
        
        if(includeTriples && !randomizeEnergies)//unrandomized search space is pairwise precomp, no triples
            throw new RuntimeException("ERROR: Can't make unrandomized search space with triples");
        
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
        
        SearchProblem ans = new SearchProblem( "examples/1CC8/testResults/CONFSEARCHTEST"+numPos, "examples/1CC8/1CC8.ss.pdb", 
                flexRes, allowedAAs,false, false, false, null, 
                false, new LUTESettings(), new DEEPerSettings(), new ArrayList<>(), new ArrayList<>(), 
                useEllipses, false, false, false, null, false, new ArrayList<>());
                //don't add WT, and no minimization, EPIC, tuple expansion, DEEPer, or strand motions

        
        if(randomizeEnergies){
            //we don't need real energies, just make some up (in fact the randomization will be good)
            ans.emat = new EnergyMatrix(ans.confSpace,0);
            ans.emat.fill(new InfiniteIterator<Double>() {
				@Override
				public Double next() {
					return getRandomEnergy();
				}
            });
            
            if(includeTriples){
                int numTriples = 5*numPos;
                
                for(int t=0; t<numTriples; t++){
                    //draw a random triple of RCs and set it to a random value
                    //we can throw in a few quadruples too actually
                    RCTuple randomTuple = randomTriple(ans,true);
                    ans.emat.setTupleValue(randomTuple, getRandomEnergy());
                }
            }
        }
        else {
            ans.loadEnergyMatrix();
        }
        
        return ans;
    }
    
    
    private static RCTuple randomTriple(SearchProblem sp, boolean allowQuad){
        //draw a random RC triple from sp
        //consider quadruples too if indicated
        Random rand = new Random();
        
        if(sp.confSpace.numPos<3)
            throw new RuntimeException("ERROR: Can't make triples in a <3-residue conf space");
        if(sp.confSpace.numPos<4)//quads impossible
            allowQuad = false;
        
        int numPos = 3;
        if(allowQuad){
            if(Math.random()>.5)//equal changes of triple or quad if allowQuad
                numPos = 4;
        }
        
        ArrayList<Integer> posList = new ArrayList<>();
        ArrayList<Integer> RCList = new ArrayList<>();
        
        for(int posCount=0; posCount<numPos; posCount++){
            
            //randomly draw position
            int pos;
            boolean posOK;
            
            do {
                posOK = true;
                pos = rand.nextInt(sp.confSpace.numPos);
                for(int count2=0; count2<posCount; count2++){//posList cannot repeat
                    if(pos == posList.get(count2))
                        posOK = false;
                }
            } while(!posOK);
            
            posList.add(pos);
            
            int rc = rand.nextInt( sp.confSpace.posFlex.get(pos).RCs.size() );
            RCList.add(rc);
        }
        
        return new RCTuple(posList, RCList);
    }
    
    private static double getRandomEnergy() {
        //randomly generate an "energy".  Keep them within a reasonable range
    	return 20*(Math.random() - 0.5);
    }
}
