/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import java.util.HashMap;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SearchProblem;

/**
 *
 * @author lowegard
 */

//takes in a list of conformations with their information and returns a dictionary 
//that contains each sequence along with it's partition function

public class PartitionFuncDict {
	
	private HashMap<String, PartitionFunction> seqDict; 
	private SearchProblem sp;
	
	public PartitionFuncDict(List<EnergiedConf> energyConfs, SearchProblem sp) {
		seqDict = new HashMap<>();
		this.sp = sp;
		for (EnergiedConf conf:energyConfs) {
			String sequence = confToSequence(conf);
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
	
	
	
	
    
}
