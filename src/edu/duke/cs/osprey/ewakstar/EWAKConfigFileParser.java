/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.StringParsing;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

/**
 *
 * @author Anna Lowegard(anna.lowegard@duke.edu)
 * Adegoke Ojewole (ao68@duke.edu)
 */
public class EWAKConfigFileParser {

	ConfigFileParser cfp;	// config file parser
	ParamSet params;		// params object
	SearchProblem[] search; // search problems

	public EWAKConfigFileParser(ConfigFileParser other) {
		this.cfp = other;
		this.params = cfp.getParams();
		this.search = null;
	}

	/**
	 * Get search problems
	 * @return
	 */
	public SearchProblem[] getSearchProblems() {
		if(search == null) {
			search = makeSearchProblems();
		}
		return search;
	}

	public void loadEnergyMatrices() {
		for(SearchProblem search : getSearchProblems()) {
			search.loadEnergyMatrix();
		}
	}

	public void pruneMatrices() {
		for(SearchProblem search : getSearchProblems()) {
			//prune
			pruneMatrix(search);
		}
	}
	
	public void pruneMatrix(SearchProblem search) {
		if(!params.getBool("UsePoissonBoltzmann")) {
			PruningControl pc = cfp.setupPruning(search, 
					params.getDouble("Ival")+params.getDouble("Ew"), 
					params.getBool("UseEpic"), 
					params.getBool("UseTupExp"));
			//silence output
			pc.setReportMode(null);
			pc.prune();
		}
	}

	/**
	 * Make and prune all search problems
	 * @return
	 */
	private SearchProblem[] makeSearchProblems() {
		int numStrands = params.getInt("NUMOFSTRANDS");
		SearchProblem[] ans = new SearchProblem[numStrands + 1];

		for (int strand = 0; strand < numStrands+1; ++strand) {
			ans[strand] = makeSearchProblem(strand, null, null);
		}

		return ans;
	}

	/**
	 * Make search problem for a specific strand
	 * @param strand
	 * @return
	 */
	public SearchProblem makeSearchProblem(int strand, ArrayList<String> mutRes, ArrayList<ArrayList<String>> allowedAAs) {
		boolean cont = params.getBool("DOMINIMIZE");
		String flexibility = cont ? "cont" : "disc";

		boolean addWT = mutRes == null ? params.getBool("AddWT") : false;
		
		if(mutRes == null) {
			mutRes = new ArrayList<>();

			int numStrands = params.searchParams("STRANDMUT").size() 
					- params.searchParams("STRANDMUTNUMS").size();

			if(strand < numStrands) {
				StringTokenizer st = new StringTokenizer(params.getValue("STRANDMUT" + strand));
				while (st.hasMoreTokens()) {
					mutRes.add(st.nextToken());
				}
			}

			else {
				mutRes = cfp.getFlexRes();
				Collections.sort(mutRes);
			}
		}

		DEEPerSettings deeperSettings = setupDEEPer(strand, mutRes, cont);
		if(!cont) deeperSettings = deeperSettings.makeDiscreteVersion();

		ArrayList<String[]> moveableStrandTermini = moveableStrandTermini(strand);
		ArrayList<String[]> freeBBZones = freeBBZoneTermini(strand);

		EPICSettings epicSettings = cont ? new EPICSettings(params) : new EPICSettings();
		LUTESettings luteSettings = cont ? new LUTESettings(params) : new LUTESettings();

		String dir = params.getValue("RunName") + File.separator + params.getValue("EmatDir");
		ObjectIO.makeDir(dir, false);

		if(allowedAAs == null) {
			allowedAAs = getAllowedAAs(mutRes);
		}
		
		SearchProblem ans = new SearchProblem(
				dir + File.separator + "Strand." + strand + "." + flexibility,
				params.getValue("PDBNAME"),
				mutRes, allowedAAs,
				addWT,
				cont,
				params.getBool("UseEPIC"),
				epicSettings,
				params.getBool("UseTupExp"),
				luteSettings,
				deeperSettings, moveableStrandTermini, freeBBZones,
				params.getBool("useEllipses"),
				params.getBool("useERef"),
				params.getBool("AddResEntropy"),
				params.getBool("addWTRots"),
				strand2Termini(strand),
				params.getBool("useVoxelG"),
				getWtRotOnlyRes()
				);

		ans.numEmatThreads = params.getInt("EmatThreads");

		return ans;
	}

	protected ArrayList<String> getWtRotOnlyRes() {
		//List of residues for which we'll only include the wild-type rotamer
		ArrayList<String> wtRotOnlyRes = new ArrayList<>();
		String val = params.getValue("WTRotOnlyRes");

		StringTokenizer tokenizer = new StringTokenizer(val);
		while (tokenizer.hasMoreTokens()) {
			wtRotOnlyRes.add(tokenizer.nextToken());
		}

		return wtRotOnlyRes;
	}

	public ArrayList<ArrayList<String>> getAllowedAAs(ArrayList<String> mutRes) {
		ArrayList<ArrayList<String>> ans = new ArrayList<>();
		for (int pos = 0; pos < mutRes.size(); ++pos) {
			ans.add(new ArrayList<>());
			String res = mutRes.get(pos);
			StringTokenizer st = new StringTokenizer(params.getValue("RESALLOWED" + res).trim());
			while (st.hasMoreTokens()) {
				ans.get(pos).add(st.nextToken());
			}
		}
		return ans;
	}

	private ResidueTermini strand2Termini(int strand) {
		ArrayList<Integer> alTmni = new ArrayList<>();
		if (params.getValue("STRAND" + strand, "").length() == 0) {
			//complex
			for (int unbound = 0; unbound < strand; ++unbound) {
				StringTokenizer st = new StringTokenizer(params.getValue("STRAND" + unbound));
				while (st.hasMoreTokens()) {
					alTmni.add(Integer.valueOf(st.nextToken()));
				}
			}
		} else {
			//unbound state
			StringTokenizer st = new StringTokenizer(params.getValue("STRAND" + strand));
			while (st.hasMoreTokens()) {
				alTmni.add(Integer.valueOf(st.nextToken()));
			}
		}
		Collections.sort(alTmni);
		return new ResidueTermini(strand, alTmni.get(0), alTmni.get(alTmni.size() - 1));
	}

	protected ArrayList<String[]> moveableStrandTermini(int subState) {
		//Read the strands that are going to translate and rotate
		//Let's say they can do this regardless of what doMinimize says (that's for sidechains)
		String key = "STRANDROTTRANS";
		ArrayList<String[]> ans = new ArrayList<>();

		for (String rt : params.searchParams(key)) {
			if (params.getBool(rt)) {
				//So rt = STRANDROTTRANS0 here means strand 0 should translate & rotate
				//OK to go through these params in lexical ordering
				String strand = rt.replaceAll(key, "").trim();

				ResidueTermini strandTmni = strand2Termini(Integer.valueOf(strand));
				ResidueTermini subStateTmni = strand2Termini(subState);

				if (!subStateTmni.contains(strandTmni)) {
					continue;
				}

				ans.add(strandTmni.toStringArray());
			}
		}
		return ans;
	}

	protected ArrayList<String[]> freeBBZoneTermini(int subState) {
		//Read the termini of the BBFreeBlocks, if any
		ArrayList<String[]> ans = new ArrayList<>();

		for (String rt : params.searchParams("BBFREEBLOCK")) {
			//So for example BBFREEBLOCK0 120 125 would mean make a BBFreeBlock for res 120-125
			//lexical ordering for blocks is OK
			String val = params.getValue(rt);

			String[] bbfTmni = {
					StringParsing.getToken(val, 1),
					StringParsing.getToken(val, 2)};

			int begin = Integer.valueOf(bbfTmni[0]);
			int end = Integer.valueOf(bbfTmni[1]);

			ResidueTermini tmni = strand2Termini(subState);
			if (tmni.contains(begin) && tmni.contains(end)) {
				ans.add(bbfTmni);
			}
		}
		return ans;
	}

	protected DEEPerSettings setupDEEPer(int strand, ArrayList<String> mutRes, boolean cont) {

		String flexibility = cont ? "cont" : "disc";

		DEEPerSettings dset = new DEEPerSettings(
				params.getBool("doPerturbations"),
				"Strand." + strand + "." + flexibility + "." + params.getRunSpecificFileName("perturbationFile", ".pert"),
				params.getBool("selectPerturbations"),
				params.getValue("startingPerturbationFile"),
				params.getBool("onlyStartingPerturbations"),
				params.getDouble("maxShearParam"),
				params.getDouble("maxBackrubParam"),
				params.getBool("selectLCAs"),
				mutRes,
				params.getValue("PDBNAME"),
				params.getBool("DORAMACHECK")
				);

		dset.loadPertFile(strand2Termini(strand));
		return dset;
	}

}
