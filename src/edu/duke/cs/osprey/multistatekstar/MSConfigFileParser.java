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

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.StringParsing;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class MSConfigFileParser extends ConfigFileParser {

	public MSConfigFileParser() {
		super();
	}

	public MSConfigFileParser(String[] confPaths) {
		super();
		for (String path : confPaths) {
			params.addParamsFromFile(path);
		}
	}

	public MSConfigFileParser(String[] confPaths, boolean isVerbose) {
		super();
		for (String path : confPaths) {
			params.addParamsFromFile(path);
		}
		params.setVerbosity(isVerbose);
	}

	public SearchProblem getSearchProblem(int state, int subState, ArrayList<Integer> mutRes, boolean cont) {

		String flexibility = cont ? "cont" : "disc";
		ArrayList<String> mutResS = new ArrayList<>();
		for(int res : mutRes) mutResS.add(String.valueOf(res));
		
		DEEPerSettings deeperSettings = setupDEEPer(state, subState, mutRes, cont);
		ArrayList<String[]> moveableUbStates = moveableUbStateTermini(subState);
		ArrayList<String[]> freeBBZones = freeBBZoneTermini(subState);

		EPICSettings epicSettings = new EPICSettings(params);
		LUTESettings luteSettings = new LUTESettings(params);

		if(!cont) {
			deeperSettings = deeperSettings.makeDiscreteVersion();
			freeBBZones = new ArrayList<>();
			moveableUbStates = new ArrayList<>();

			epicSettings = new EPICSettings();
			luteSettings = new LUTESettings();
		}

		String dir = params.getValue("RunName")+File.separator+params.getValue("EmatDir");
		ObjectIO.makeDir(dir, false);
		
		SearchProblem ans = new SearchProblem(
				dir+File.separator+"State."+state+"."+subState+"."+flexibility, 
				params.getValue("PDBNAME"), 
				mutResS, getAllowedAAs(mutRes),
				params.getBool("AddWT"), 
				cont,
				params.getBool("UseEPIC"),
				epicSettings,
				params.getBool("UseTupExp"),
				luteSettings,
				deeperSettings, moveableUbStates, freeBBZones,
				params.getBool("useEllipses"),
				params.getBool("useERef"),
				params.getBool("AddResEntropy"),
				params.getBool("addWTRots"),
				subState2Termini(subState),
				params.getBool("useVoxelG"),
				getWtRotOnlyRes()
				);

		ans.numEmatThreads = params.getInt("EmatThreads");

		return ans;
	}

	public String[] getWtSeq(ArrayList<Integer> mutRes) {
		Molecule m = new Strand.Builder(PDBIO.readFile(params.getValue("PDBNAME"))).build().mol;
		int numPos = mutRes.size();
		String[] wt = new String[numPos];

		for(int pos=0; pos<numPos; pos++) {
			Residue res = m.getResByPDBResNumber(String.valueOf(mutRes.get(pos)));
			String wtName = res.template.name;
			wt[pos] = wtName;
		}
		return wt;
	}

	public ArrayList<ArrayList<String>> getAllowedAAs(ArrayList<Integer> mutRes) {
		ArrayList<ArrayList<String>> ans = new ArrayList<>();
		for(int pos=0;pos<mutRes.size();++pos) {
			ans.add(new ArrayList<>());
			int res = mutRes.get(pos);
			StringTokenizer st = new StringTokenizer(params.getValue("RESALLOWED"+res).trim());
			while(st.hasMoreTokens()) ans.get(pos).add(st.nextToken());
		}
		return ans;
	}

	ResidueTermini subState2Termini(int subState) {
		String key = params.searchParams("UbStateLimits").size() > 0 ? "UbStateLimits" : "Strand";
		ArrayList<String> alTmni = new ArrayList<>();
		if(params.getValue(key+subState, "").length()==0) {
			//complex
			for(int unbound=0;unbound<subState;++unbound) {
				StringTokenizer st = new StringTokenizer(params.getValue(key+unbound));
				while(st.hasMoreTokens()) alTmni.add(st.nextToken());
			}
		}
		else {
			//unbound state
			StringTokenizer st = new StringTokenizer(params.getValue(key+subState));
			while(st.hasMoreTokens()) alTmni.add(st.nextToken());
		}
		Collections.sort(alTmni);
		return new ResidueTermini(subState, alTmni.get(0), alTmni.get(alTmni.size()-1));
	}

	protected DEEPerSettings setupDEEPer(int state, int subState, ArrayList<Integer> mutRes, boolean cont) {

		ArrayList<String> sMutRes = new ArrayList<>();
		for(int res : mutRes) sMutRes.add(String.valueOf(res));
		String flexibility = cont ? "cont" : "disc";

		DEEPerSettings dset = new DEEPerSettings(
				params.getBool("doPerturbations"),
				"State."+state+"."+subState+"."+flexibility+"."+params.getRunSpecificFileName("perturbationFile", ".pert"),
				params.getBool("selectPerturbations"),
				params.getValue("startingPerturbationFile"),
				params.getBool("onlyStartingPerturbations"),
				params.getDouble("maxShearParam"),
				params.getDouble("maxBackrubParam"),
				params.getBool("selectLCAs"),
				sMutRes,
				params.getValue("PDBNAME"),
				params.getBool("DORAMACHECK"),
				EnvironmentVars.resTemplates
				);

		dset.loadPertFile(subState2Termini(subState));
		return dset;
	}

	protected ArrayList<String[]> freeBBZoneTermini(int subState) {
		//Read the termini of the BBFreeBlocks, if any
		ArrayList<String[]> ans = new ArrayList<>();

		for(String rt : params.searchParams("BBFREEBLOCK")){
			//So for example BBFREEBLOCK0 120 125 would mean make a BBFreeBlock for res 120-125
			//lexical ordering for blocks is OK
			String val = params.getValue(rt);

			String[] bbfTmni = { 
					StringParsing.getToken(val, 1), 
					StringParsing.getToken(val, 2) };

			ResidueTermini tmni = subState2Termini(subState);
			if(tmni.contains(bbfTmni[0]) && tmni.contains(bbfTmni[1]))
				ans.add(bbfTmni);
		}
		return ans;
	}

	protected ArrayList<String[]> moveableUbStateTermini(int subState) {
		//Read the strands that are going to translate and rotate
		//Let's say they can do this regardless of what doMinimize says (that's for sidechains)
		String key = params.searchParams("UBSTATEROTTRANS").size() > 0 ? "UBSTATEROTTRANS" : "STRANDROTTRANS";
		ArrayList<String[]> ans = new ArrayList<>();
		for(String rt : params.searchParams(key)){
			if(params.getBool(rt)){
				//So rt = UBSTATEROTTRANS0 here means strand 0 should translate & rotate
				//OK to go through these params in lexical ordering
				String ubState = rt.replaceAll(key, "").trim();
				if(!String.valueOf(subState).equals(ubState)) continue;
				ResidueTermini tmni = subState2Termini(subState);
				ans.add(tmni.toStringArray());
			}
		}
		return ans;
	}

}
