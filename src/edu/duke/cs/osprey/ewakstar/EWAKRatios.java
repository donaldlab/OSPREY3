/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.GMECFinder;

/**
 *
 * @author Anna Lowegard(anna.lowegard@duke.edu)
 * Adegoke Ojewole (ao68@duke.edu)
 */
public class EWAKRatios {

	private ConfigFileParser cfp;
	private PartitionFuncDict complexPfd;
	private ArrayList<PartitionFuncDict> strandPfd;
	private HashMap<Integer, HashSet<String>> unboundAllowedSeqsByStrand;
	private HashMap<Integer, HashSet<String>> missedSeqsByStrand;
	private int numStrands;
	
	public EWAKRatios(ConfigFileParser cfp) {
		
		this.cfp = cfp;
		this.strandPfd = null;
		this.unboundAllowedSeqsByStrand = null;
		this.missedSeqsByStrand = null;
		
		this.numStrands = cfp.getParams().searchParams("STRANDMUT").size() 
				- cfp.getParams().searchParams("STRANDMUTNUMS").size();
		
		GMECFinder complexes = new GMECFinder();
		complexes.init(cfp);
		List<EnergiedConf> complexConfs = complexes.calcGMEC();

		complexPfd = new PartitionFuncDict(complexConfs, 
				cfp.getSearchProblem(),
				null);
		
		//get allowed unbound sequences
		unboundAllowedSeqsByStrand = getUnboundAllowedSeqsByStrand();
		//make unbound state partition function dictionaries
		strandPfd = createUnboundPartitionFuncDicts();
		missedSeqsByStrand = computeMissedSeqsByStrand();
		
		System.out.println();
		System.out.println("Printing complex sequences ...");
		complexPfd.printSequences();
		System.out.println(" ... done!");
		
		for(int strand = 0; strand < strandPfd.size(); ++strand) {
			PartitionFuncDict pfd = strandPfd.get(strand);
			System.out.println();
			System.out.println("Printing strand"+strand+" sequences ...");
			pfd.printSequences();
			//print missed sequences
			System.out.println("Missed sequences");
			for(String seq : missedSeqsByStrand.get(strand)) {
				System.out.println(seq);
			}
			System.out.println(" ... done!");
		}

		/* TODO:
		 * 1) iterate through list of complex confs
		 * 		map confs to sequence
		 * 		update partition function of sequence
		 * 2) make list of P and L only sequences from complex sequences
		 * 		limit P and L pruning matrices accordingly
		 * 3) repeat step 1 for P and L
		 * 4) print output to file
		 */


		/*
    	// parse config file
    	EWAKConfigFileParser ecfp = new EWAKConfigFileParser(cfp);
    	// make search problem
        SearchProblem[] sps = ecfp.getSearchProblems();
        ecfp.loadEnergyMatrices();
        ecfp.pruneMatrices();

        GMECFinder[] gfs = new GMECFinder[sps.length];
        for(int i = 0; i < gfs.length; ++i) {
        	gfs[i] = new GMECFinder();
        	gfs[i].init(cfp, sps[i]);
        	gfs[i].calcGMEC();
        }
		*/
	}

	public ArrayList<PartitionFuncDict> createUnboundPartitionFuncDicts() {
		ArrayList<PartitionFuncDict> ans = new ArrayList<>();
		EWAKConfigFileParser ecfp = new EWAKConfigFileParser(cfp);

		//creating strand searchproblems+emats
		for(int strand = 0; strand < numStrands; ++strand) {
			//get res2allowedaas
			LinkedHashMap<Integer, Set<String>> res2AllowedAAs = complexPfd.getAllowedAAs(cfp, strand);

			ArrayList<String> mutRes = new ArrayList<>();
			ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();

			//convert allowed aas to form accepted by config file parser
			for(Integer res : res2AllowedAAs.keySet()) {
				
				System.out.print("RESALLOWED"+res+" ");
				StringBuilder sb = new StringBuilder();
				for(String aa : res2AllowedAAs.get(res)) {
					sb.append(aa + " ");
				}
				System.out.println(sb.toString().trim());
				
				mutRes.add(res.toString());
				allowedAAs.add(new ArrayList<>(res2AllowedAAs.get(res)));

			}
			
			SearchProblem search = ecfp.makeSearchProblem(strand, mutRes, allowedAAs);
			search.loadEnergyMatrix();
			ecfp.pruneMatrix(search);
			
			GMECFinder strandGMEC = new GMECFinder();
			strandGMEC.init(cfp, search);
			List<EnergiedConf> strandConfs = strandGMEC.calcGMEC();
			ans.add(new PartitionFuncDict(strandConfs, 
					search, 
					unboundAllowedSeqsByStrand.get(strand)));
		}
		
		return ans;
	}
	
	private HashMap<Integer, HashSet<String>> getUnboundAllowedSeqsByStrand() {
		HashMap<Integer, HashSet<String>> ans = new HashMap<>();
		
		for(int strand = 0; strand < numStrands; ++strand) {
			SubSequenceBounds ssb = new SubSequenceBounds(cfp, strand);
			
			HashSet<String> strandSeqs = new HashSet<>();
			for(String seq : complexPfd.getSequences()) {
				String strandSeq = seq.substring(ssb.start, Math.min(ssb.end, seq.length())).trim();
				strandSeqs.add(strandSeq);
			}
			
			ans.put(strand, strandSeqs);
		}
		
		return ans;
	}
	
	private HashMap<Integer, HashSet<String>> computeMissedSeqsByStrand() {
		HashMap<Integer, HashSet<String>> ans = new HashMap<>();
		for(int strand = 0; strand < numStrands; ++strand) {
			ans.put(strand, new HashSet<>());
			for(String seq : unboundAllowedSeqsByStrand.get(strand)) {
				if(!strandPfd.get(strand).getSequences().contains(seq)) {
					ans.get(strand).add(seq);
				}
			}
		}
		return ans;
	}
	
	public void run() {
	};

}
