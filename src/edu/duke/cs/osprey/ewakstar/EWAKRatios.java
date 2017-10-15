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
import java.util.StringTokenizer;

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
		//we need type dependent DEE to be true
		cfp.getParams().setValue("TYPEDEP", "TRUE");

		this.cfp = cfp;
		this.strandPfd = null;
		this.unboundAllowedSeqsByStrand = null;
		this.missedSeqsByStrand = null;

		this.numStrands = cfp.getParams().searchParams("STRANDMUT").size() 
				- cfp.getParams().searchParams("STRANDMUTNUMS").size();
	}

	public ArrayList<PartitionFuncDict> createUnboundPartitionFuncDicts() {
		ArrayList<PartitionFuncDict> ans = new ArrayList<>();

		//creating strand searchproblems+emats
		for(int strand = 0; strand < numStrands; ++strand) {
			//get res2allowedaas
			LinkedHashMap<Integer, Set<String>> res2AllowedAAs = complexPfd.getAllowedAAs(cfp, strand);

			ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();

			//convert allowed aas to form accepted by config file parser
			for(Integer res : res2AllowedAAs.keySet()) {

				System.out.print("RESALLOWED"+res+" ");
				StringBuilder sb = new StringBuilder();
				for(String aa : res2AllowedAAs.get(res)) {
					sb.append(aa + " ");
				}
				System.out.println(sb.toString().trim());

				allowedAAs.add(new ArrayList<>(res2AllowedAAs.get(res)));
			}

			ans.add(createUnboundPartitionFuncDict(strand, allowedAAs, unboundAllowedSeqsByStrand.get(strand)));
		}

		return ans;
	}

	private PartitionFuncDict createUnboundPartitionFuncDict(int strand, 
			ArrayList<ArrayList<String>> allowedAAs,
			HashSet<String> allowedSeqs) {

		EWAKConfigFileParser ecfp = new EWAKConfigFileParser(cfp);
		SearchProblem search = ecfp.makeSearchProblem(strand);
		search.loadEnergyMatrix();
		ecfp.pruneMatrix(search);
		
		ArrayList<Integer> pos = new ArrayList<>();
		for(int i = 0; i < allowedAAs.size(); ++i) pos.add(i);
		EWAKSearchProblem ewakSearch = new EWAKSearchProblem(search, pos, allowedAAs);
		ewakSearch.updatePruningMatrix();

		GMECFinder strandGMEC = new GMECFinder();
		strandGMEC.init(cfp, ewakSearch);
		List<EnergiedConf> strandConfs = strandGMEC.calcGMEC();
		return new PartitionFuncDict(strandConfs, ewakSearch, allowedSeqs);
	}

	public HashMap<Integer, ArrayList<PartitionFuncDict>> computeMissedUnboundPartitionFuncDicts() {
		HashMap<Integer, ArrayList<PartitionFuncDict>> ans = new HashMap<>();

		for(int strand = 0; strand < numStrands; ++strand) {
			ans.put(strand, new ArrayList<>());

			for(String strandSeq : missedSeqsByStrand.get(strand)) {
				ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();

				StringTokenizer st = new StringTokenizer(strandSeq);
				while(st.hasMoreTokens()) {
					ArrayList<String> aaAtPos = new ArrayList<>();
					aaAtPos.add(st.nextToken().split("-")[0].trim());
					allowedAAs.add(aaAtPos);
				}

				HashSet<String> singleSeq = new HashSet<>();
				singleSeq.add(strandSeq);

				ans.get(strand).add(createUnboundPartitionFuncDict(strand, allowedAAs, singleSeq));
			}
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

	private HashMap<Integer, HashSet<String>> getMissedSeqsByStrand() {
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

	void printMetaData() {
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
			if(missedSeqsByStrand.containsKey(strand)) {
				for(String seq : missedSeqsByStrand.get(strand)) {
					System.out.println(seq);
				}
			}
			System.out.println(" ... done!");
		}
	}

	public void run() {
		/* ALGORITHM:
		 * 1) iterate through list of complex confs
		 * 		map confs to sequence
		 * 		update partition function of sequence
		 * 2) make list of P and L only sequences from complex sequences
		 * 		limit P and L pruning matrices accordingly
		 * 3) repeat step 1 for P and L
		 * 		catch stragglers: sequences that are not enumerated by ival+ew
		 * 4) print output to file
		 */

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
		//computed sequences not enumeated by ival+ew
		missedSeqsByStrand = getMissedSeqsByStrand();

		printMetaData();

		//if ival+ew missed any unbound sequences, compute them individually and
		//merge them into the correct strand partition function dictionary
		for(HashSet<String> seqsByStrand : missedSeqsByStrand.values()) {
			if(seqsByStrand.size() > 0) {

				//compute missed sequences
				HashMap<Integer, ArrayList<PartitionFuncDict>> missedSeqDicts = computeMissedUnboundPartitionFuncDicts();
				missedSeqsByStrand.clear();

				//merge dictionaries
				for(Integer strand : missedSeqDicts.keySet()) {
					for(PartitionFuncDict dict : missedSeqDicts.get(strand)) {
						strandPfd.get(strand).merge(dict);
					}
				}

				printMetaData();

				break;
			}
		}
	}

}
