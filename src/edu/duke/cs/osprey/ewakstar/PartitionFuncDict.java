/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;

/**
 *
 * @author Anna Lowegard(anna.lowegard@duke.edu)
 * Adegoke Ojewole (ao68@duke.edu)
 */

//takes in a list of conformations with their information and returns a dictionary 
//that contains each sequence along with it's partition function

public class PartitionFuncDict {
	
	private HashMap<String, PartitionFunction> seqDict; 
	private SearchProblem sp;
	
	public PartitionFuncDict(List<EnergiedConf> energyConfs, SearchProblem sp, HashSet<String> allowedSeqs) {
		seqDict = new HashMap<>();
		this.sp = sp;
		for (EnergiedConf conf:energyConfs) {
			String sequence = confToSequence(conf);
			
			//only add allowed sequences to dictionary
			if(allowedSeqs != null && !allowedSeqs.contains(sequence)) {
				continue;
			}
			
			//if the sequence is in the dictionary, update the pf with the new energy
			if (seqDict.containsKey(sequence)) { 
				seqDict.get(sequence).addEnergy(conf.getEnergy());
			}
			
			//if the sequence is not in the dictionary, create a new partition function,
			//add the initial energy to the pf, and add the sequence & pf to the dict
			else {
				PartitionFunction pf = sequenceToPartitionFunction(sequence);
				pf.addEnergy(conf.getEnergy());
				seqDict.put(sequence, pf);
			}
			
		}
	}
	
	//import partition function dictionary entries from other
	public void merge(PartitionFuncDict other) {
		for(String seq : other.seqDict.keySet()) {
			if(this.seqDict.containsKey(seq)) {
				System.out.println("WARNING: sequence " + seq + " already exists. skipping...");
				continue;
			}
			
			this.seqDict.put(seq, other.seqDict.get(seq));
		}
	}
	
	//constructs a new partition function for a given sequence
	public PartitionFunction sequenceToPartitionFunction(String sequence) {
		PartitionFunction pf = new PartitionFunction();
		pf.setSequence(sequence);
		return pf;
	}
	
	//makes a user readable version of the sequence e.g. ARG-39
	public String confToSequence(EnergiedConf conf) {
		StringBuilder sb = new StringBuilder();
		int length = conf.getAssignments().length;
		for (int i=0; i<length; i++) {
			int rotAssignment = conf.getAssignments()[i];
			String aminoAcid = sp.confSpace.posFlex.get(i).RCs.get(rotAssignment).AAType;
			String residue = sp.flexRes.get(i);
			sb.append(aminoAcid+"-"+residue+" ");
		}
		
		String ans = sb.toString().trim();
		return ans;
	}
	
	//gets sequences in dictionary
	public Set<String> getSequences() {
		return seqDict.keySet();
	}
	
	//print formatted sequences
	public void printSequences() {
		for(String key : seqDict.keySet()) {
			System.out.println(seqDict.get(key).toString());
		}
	}
	
	//filters allowed aas by residue
	public LinkedHashMap<Integer, Set<String>> getAllowedAAsByResidue() {
		LinkedHashMap<Integer, Set<String>> ans = new LinkedHashMap<>();
		
		for(String seq : seqDict.keySet()) {
			StringTokenizer st = new StringTokenizer(seq);
			while(st.hasMoreTokens()) {
				String token = st.nextToken();
				int res = Integer.valueOf(token.split("-")[1]);
				String aa = token.split("-")[0];
				
				if(!ans.containsKey(res)) {
					ans.put(res, new HashSet<String>());
				}
				
				ans.get(res).add(aa);
			}
		}
		
		return ans;
	}
	
	//filters allowed aas by strand
	public LinkedHashMap<Integer, Set<String>> getAllowedAAs(ConfigFileParser cfp, int strand) {
		//get set of flexible residues in specified strand
		String flexResByStrand = cfp.getParams().getValue("STRANDMUT"+strand, "").trim();
		HashSet<Integer> flexResSet = new HashSet<>();
		StringTokenizer st = new StringTokenizer(flexResByStrand);
		while(st.hasMoreTokens()) {
			flexResSet.add(Integer.valueOf(st.nextToken()));
		}
		
		//remove residues not in this set
		LinkedHashMap<Integer, Set<String>> all = getAllowedAAsByResidue();
		HashSet<Integer> toRemove = new HashSet<>();
		for(Integer res : all.keySet()) {
			if(!flexResSet.contains(res)) {
				toRemove.add(res);
			}
		}
		
		for(Integer key : toRemove) all.remove(key);
		return all;
	}
    
}
