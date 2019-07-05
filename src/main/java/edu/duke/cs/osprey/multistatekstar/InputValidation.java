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

package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.restypes.DAminoAcidHandler;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.StringParsing;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */
public class InputValidation {

	ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions;
	ArrayList<ArrayList<ArrayList<Integer>>> state2MutableResNums;
	
	public InputValidation(ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions,
			ArrayList<ArrayList<ArrayList<Integer>>> state2MutableResNums) {
		this.AATypeOptions = AATypeOptions;
		this.state2MutableResNums = state2MutableResNums;
	}
	
	public void handleObjFcn(ParamSet msParams, LMB objFcn) {
		if(objFcn.getCoeffs().length != msParams.getInt("NUMSTATES"))
			throw new RuntimeException("ERROR: the number of OBJFCN coefficients must equal NUMSTATES");
		for(BigDecimal coeff : objFcn.getCoeffs())
			if(coeff.compareTo(BigDecimal.ZERO)==0)
				throw new RuntimeException("ERROR: objective function coefficient cannot be 0");
	}
	
	public void handleConstraints(ParamSet msParams, LMB[] constraints) {
		for(LMB constr : constraints) {
			if(constr.getCoeffs().length != msParams.getInt("NUMSTATES"))
				throw new RuntimeException("ERROR: the number of constraint coefficients must equal NUMSTATES");
		}
	}
	
	public void handleAATypeOptions(int state, int subState, MSConfigFileParser stateCfp) {
		//Given the config file parser for a state, make sure AATypeOptions
		//matches the allowed AA types for this state

		Molecule wtMolec = new Strand.Builder(PDBIO.readFile(stateCfp.params.getValue("PDBName"))).build().mol;

		ArrayList<Integer> mutRes = state2MutableResNums.get(state).get(subState);
		ArrayList<ArrayList<String>> subStateAAOptions = stateCfp.getAllowedAAs(mutRes);

		if(AATypeOptions==null) 
			AATypeOptions = new ArrayList<>();

		for(int mutPos=0; mutPos<mutRes.size(); mutPos++) {
			ArrayList<String> subStateResOptions = subStateAAOptions.get(mutPos);
			if(stateCfp.params.getBool("AddWT")) {
				Residue res = wtMolec.getResByPDBResNumber(String.valueOf(mutRes.get(mutPos)));
				//always add wt in pos 0
				if(StringParsing.containsIgnoreCase(subStateResOptions, res.template.name))
					subStateResOptions.remove(res.template.name);
				subStateResOptions.add(0, res.template.name);
			}

			if(AATypeOptions.size()<=state)//add storage for state
				AATypeOptions.add(new ArrayList<>());

			if(AATypeOptions.get(state).size()<=subState)//add storage for substate
				AATypeOptions.get(state).add(new ArrayList<>());

			if(AATypeOptions.get(state).get(subState).size()<=mutPos)//need to fill in based on this substate
				AATypeOptions.get(state).get(subState).add(subStateResOptions);

			//check for correspondence in AAs among states
			ArrayList<String> defaultSubStateResOptions = AATypeOptions.get(0).get(subState).get(mutPos);

			if(defaultSubStateResOptions.size()!=subStateResOptions.size()){
				throw new RuntimeException("ERROR: Current state has "+
						subStateResOptions.size()+" AA types allowed for position "+mutPos
						+" compared to "+defaultSubStateResOptions.size()+" for previous states");
			}

			for(int a=0;a<defaultSubStateResOptions.size();++a){
				String aa1 = defaultSubStateResOptions.get(a);
				String aa2 = subStateResOptions.get(a);

				//only amino acids must correspond between states
				if(!DAminoAcidHandler.isStandardLAminoAcid(aa1) || 
						!DAminoAcidHandler.isStandardLAminoAcid(aa2)) continue;

				if(!aa1.equalsIgnoreCase(aa2))
					throw new RuntimeException("ERROR: Current state has AA type "+
							aa2+" where previous states have AA type "+aa1+", at position "+mutPos);
			}
		}
	}
	
	public void handleStateParams(int state, ParamSet sParams, ParamSet msParams) {
		//parameter sanity check
		boolean imindee = sParams.getBool("IMINDEE");
		boolean doMinimize = sParams.getBool("DOMINIMIZE");
		if(imindee != doMinimize)
			throw new RuntimeException("ERROR: IMINDEE must have the same value as DOMINIMIZE");
		
		double epsilon = sParams.getDouble("EPSILON");
		if(epsilon >= 1 || epsilon < 0) 
			throw new RuntimeException("ERROR: EPSILON must be >= 0 and < 1"); 
		
		//check number of constraints
		int numUbConstr = sParams.getInt("NUMUBCONSTR");
		ArrayList<String> ubConstr = sParams.searchParams("UBCONSTR");
		ubConstr.remove("NUMUBCONSTR");
		if(numUbConstr != ubConstr.size())
			throw new RuntimeException("ERROR: NUMUBCONSTR != number of listed constraints");

		int numUbStates = sParams.getInt("NUMUBSTATES");
		if(numUbStates<2) throw new RuntimeException("ERROR: NUMUBSTATES must be >=2");

		String ubStateMutNums = sParams.getValue("UBSTATEMUTNUMS");
		StringTokenizer st = new StringTokenizer(ubStateMutNums);
		if(st.countTokens() != numUbStates) throw new RuntimeException("ERROR: "
				+ "the number of tokens in UBSTATEMUTNUMS should be the same as "
				+ "NUMUBSTATES");

		if(numUbStates!=sParams.searchParams("UBSTATELIMITS").size())
			throw new RuntimeException("ERROR: need an UBSTATELIMITS line for each NUMUBSTATES");

		int numMutsRes = 0;
		while(st.hasMoreTokens()) numMutsRes += Integer.valueOf(st.nextToken());
		if(numMutsRes != msParams.getInt("NUMMUTRES")) throw new RuntimeException("ERROR: "
				+"UBSTATEMUTNUMS does not sum up to NUMMUTRES");

		ArrayList<Integer> globalMutList = new ArrayList<>();
		for(int ubState=0;ubState<numUbStates;++ubState) {
			//num unbound state residues must match number of listed residues
			int numUbMutRes = Integer.valueOf(StringParsing.getToken(ubStateMutNums, ubState+1));
			String ubMutRes = sParams.getValue("UBSTATEMUT"+ubState);
			st = new StringTokenizer(ubMutRes);
			ArrayList<Integer> ubStateMutList = new ArrayList<>();
			while(st.hasMoreTokens()) ubStateMutList.add(Integer.valueOf(st.nextToken()));
			ubStateMutList = new ArrayList<>(new HashSet<>(ubStateMutList));
			if(ubStateMutList.size()!=numUbMutRes) throw new RuntimeException("ERROR: the "
					+"number of distinct mutable residues in UBSTATEMUT"+ubState+
					" is not equal to the value specified in UBSTATEMUTNUMS");

			globalMutList.addAll(ubStateMutList);

			//listed unbound state residues must be within limits
			ArrayList<Integer> ubStateLims = new ArrayList<>();
			st = new StringTokenizer(sParams.getValue("UBSTATELIMITS"+state));
			if(st.countTokens()!=2) throw new RuntimeException("ERROR: UBSTATELIMITS"+state
					+" must have 2 tokens");
			while(st.hasMoreTokens()) ubStateLims.add(Integer.valueOf(st.nextToken()));
			Collections.sort(ubStateLims);
			for(int res : ubStateMutList){
				if(res<ubStateLims.get(0) && res>ubStateLims.get(1)) throw new RuntimeException("ERROR: "
						+"mutable residue "+res+" exceeds the boundaries of UBSTATELIMITS"+state);
			}

			//ResAllowed must exist for each mutable residue
			for(int res: ubStateMutList) {
				if(sParams.getValue("RESALLOWED"+res, "").length()==0)
					throw new RuntimeException("ERROR: RESALLOWED"+res+" must be delcared");
			}
		}

		//check that all RESALLOWED is in list of mutable residues
		ArrayList<String> raKeys = sParams.searchParams("RESALLOWED");
		if(raKeys.size() > numMutsRes) {
			for(String raVal : raKeys) {
				raVal = raVal.replaceAll("RESALLOWED", "").trim();
				if(!globalMutList.contains(Integer.valueOf(raVal)))
					throw new RuntimeException("ERROR: RESALLOWED"+raVal+" is not in the list of UBSTATEMUT");
			}
		}

		//check that ubState limits are mutually exclusive 
		ArrayList<ArrayList<Integer>> ubStateLimits = new ArrayList<>();
		for(int ubState=0;ubState<numUbStates;++ubState) {
			st = new StringTokenizer(sParams.getValue("UBSTATELIMITS"+ubState));
			ArrayList<Integer> tmp = new ArrayList<Integer>();
			while(st.hasMoreTokens()) tmp.add(Integer.valueOf(st.nextToken()));
			Collections.sort(tmp);
			ubStateLimits.add(tmp);
		}
		if(ubStateLimits.get(0).get(0) <= ubStateLimits.get(1).get(1) && 
				ubStateLimits.get(1).get(0) <= ubStateLimits.get(0).get(1))
			throw new RuntimeException("ERROR: UBSTATELIMITS are not disjoint");
	}
	
}
