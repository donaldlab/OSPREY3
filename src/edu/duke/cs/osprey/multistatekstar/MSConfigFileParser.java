package edu.duke.cs.osprey.multistatekstar;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.StringParsing;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class MSConfigFileParser extends ConfigFileParser {

	public MSConfigFileParser(String[] args, boolean isVerbose) {
		super(args, isVerbose);
	}

	public MSConfigFileParser(String[] args) {
		super(args);
	}

	public MSConfigFileParser() {
		super();
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
		Molecule m = PDBFileReader.readPDBFile(params.getValue("PDBNAME"), null);
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
		String key = "STRAND";
		ArrayList<Integer> alTmni = new ArrayList<>();
		if(params.getValue(key+subState, "").length()==0) {
			//complex
			for(int unbound=0;unbound<subState;++unbound) {
				StringTokenizer st = new StringTokenizer(params.getValue(key+unbound));
				while(st.hasMoreTokens()) alTmni.add(Integer.valueOf(st.nextToken()));
			}
		}
		else {
			//unbound state
			StringTokenizer st = new StringTokenizer(params.getValue(key+subState));
			while(st.hasMoreTokens()) alTmni.add(Integer.valueOf(st.nextToken()));
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
				params.getBool("DORAMACHECK")
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

			int begin = Integer.valueOf(bbfTmni[0]);
			int end = Integer.valueOf(bbfTmni[1]);

			ResidueTermini tmni = subState2Termini(subState);
			if(tmni.contains(begin) && tmni.contains(end))
				ans.add(bbfTmni);
		}
		return ans;
	}

	protected ArrayList<String[]> moveableUbStateTermini(int subState) {
		//Read the strands that are going to translate and rotate
		//Let's say they can do this regardless of what doMinimize says (that's for sidechains)
		String key = "STRANDROTTRANS";
		ArrayList<String[]> ans = new ArrayList<>();
		for(String rt : params.searchParams(key)){
			if(params.getBool(rt)){
				//So rt = STRANDROTTRANS0 here means strand 0 should translate & rotate
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
