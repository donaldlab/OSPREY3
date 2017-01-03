package edu.duke.cs.osprey.kstar;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
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

public class KSConfigFileParser extends ConfigFileParser implements Serializable {

	private static final long serialVersionUID = 3653519769830569960L;


	public KSConfigFileParser(String[] args) {
		super(args);
	}
	
	
	public KSConfigFileParser() {
		super();
	}


	public ArrayList<ArrayList<String>> getHighOrderTuplesByStrand(int strand) {

		ArrayList<ArrayList<String>> ans = getHighOrderTuplesByPDBResNum();

		switch(strand) {

		case KSTermini.COMPLEX:
			return ans;

		case KSTermini.LIGAND:
			return filterHOTListByStrand(getStrandLimits(KSTermini.LIGAND), ans);

		case KSTermini.PROTEIN:
			return filterHOTListByStrand(getStrandLimits(KSTermini.PROTEIN), ans);

		default:
			throw new RuntimeException("ERROR: invalid strand");
		}
	}


	public ArrayList<ArrayList<String>> getHighOrderTuplesByPDBResNum() {
		ArrayList<ArrayList<String>> ans = new ArrayList<>();

		HashSet<String> flexRes = new HashSet<>(getFlexRes());
		HashSet<String> resInHot = new HashSet<>();

		ArrayList<String> hotRecords = params.searchParams("kStarPFuncHot_");

		for(String hotRecord : hotRecords) {
			ArrayList<String> hotRes = new ArrayList<>();

			String val = params.getValue(hotRecord);
			StringTokenizer tokenizer = new StringTokenizer(val);

			while(tokenizer.hasMoreTokens()) {
				String token = tokenizer.nextToken();

				if(!flexRes.contains(token))
					throw new RuntimeException("ERROR: residue " + token + " in line " + hotRecord + " must be flexible.");

				if(resInHot.contains(token))
					throw new RuntimeException("ERROR: residue " + token + " cannot appear in more than one HOT.");

				hotRes.add(token);
				resInHot.add(token);
			}

			if(hotRes.size() < 3)
				throw new RuntimeException("ERROR: in line " + hotRecord + " the number of residues is < 3.");

			hotRes.trimToSize();
			ans.add(hotRes);
		}

		ans.trimToSize();
		return ans;
	}


	protected ArrayList<ArrayList<String>> filterHOTListByStrand( KSTermini strand, ArrayList<ArrayList<String>> list ) {
		@SuppressWarnings("unchecked")
		ArrayList<ArrayList<String>> ans = (ArrayList<ArrayList<String>>) ObjectIO.deepCopy(list);

		for( Iterator<ArrayList<String>> iterator = ans.iterator(); iterator.hasNext(); ) {

			ArrayList<String> hot = iterator.next();

			for( Iterator<String> iterator2 = hot.iterator(); iterator2.hasNext(); ) {
				String pdbResNum = iterator2.next();

				if(!strand.contains(Integer.parseInt(pdbResNum)))
					iterator2.remove();
			}

			if(hot.size() < 3)
				iterator.remove();
		}

		return ans;
	}


	protected DEEPerSettings setupDEEPer(int strand) {
		//Set up the DEEPerSettings object, including the PertSet (describes the perturbations)
		//String runName = params.getValue("runName");

		DEEPerSettings dset = new DEEPerSettings(
				params.getBool("doPerturbations"),
				"STR"+strand+"."+params.getRunSpecificFileName("perturbationFile", ".pert"),
				params.getBool("selectPerturbations"),
				params.getValue("startingPerturbationFile"),
				params.getBool("onlyStartingPerturbations"),
				params.getDouble("maxShearParam"),
				params.getDouble("maxBackrubParam"),
				params.getBool("selectLCAs"),
				getFlexResByStrand(strand),
				params.getValue("PDBNAME"),
				params.getBool("DORAMACHECK")
				);

		// remove residues not in this strand
		KSTermini limits = getStrandLimits(strand);

		// perturbation file is by strand
		ObjectIO.delete(params.getRunSpecificFileName("perturbationFile", ".pert"));

		dset.loadPertFile(limits);//load the PertSet from its file
		return dset;
	}


	private ArrayList<String[]> freeBBZoneTermini(KSTermini limits){
		//Read the termini of the BBFreeBlocks, if any
		ArrayList<String[]> ans = new ArrayList<>();

		for(String rt : params.searchParams("BBFREEBLOCK")){
			//So for example BBFREEBLOCK0 120 125 would mean make a BBFreeBlock for res 120-125
			//lexical ordering for blocks is OK
			String strandLimitsString = params.getValue(rt);

			String[] termini = 
				{ StringParsing.getToken(strandLimitsString, 1),
						StringParsing.getToken(strandLimitsString, 2) };

			int begin = Integer.parseInt(termini[0]);
			int end = Integer.parseInt(termini[1]);

			if(limits == null || (limits.contains(begin) && limits.contains(end)))
				ans.add(termini);
		}

		return ans;
	}


	private ArrayList<String[]> moveableStrandTermini(KSTermini limits){
		//Read the strands that are going to translate and rotate
		//Let's say they can do this regardless of what doMinimize says (that's for sidechains)
		ArrayList<String[]> ans = new ArrayList<>();

		for(String rt : params.searchParams("STRANDROTTRANS")){
			if(params.getBool(rt)){
				//So rt = STRANDROTTRANS0 here means strand 0 should translate & rotate
				//OK to go through these params in lexical ordering
				String strandNum = rt.substring(14);
				String strandLimitsString = params.getValue("STRAND"+strandNum);

				String[] termini = 
					{ StringParsing.getToken(strandLimitsString, 1),
							StringParsing.getToken(strandLimitsString, 2) };

				int begin = Integer.parseInt(termini[0]);
				int end = Integer.parseInt(termini[1]);

				if(limits == null || (limits.contains(begin) && limits.contains(end)))
					ans.add(termini);
			}
		}

		return ans;
	}


	protected KSTermini getStrandLimits( int strand ) {

		if( strand == KSTermini.COMPLEX ) return null;

		String strandLimits = params.getValue( "STRAND"+strand );

		StringTokenizer tokenizer = 
				new StringTokenizer(params.getValue( "STRANDMUTNUMS" ));
		for( int it = 0; it < strand; ++it ) tokenizer.nextToken(); 
		int numFlexRes = Integer.parseInt( tokenizer.nextToken() );

		tokenizer = new StringTokenizer( strandLimits );

		ArrayList<String> strLimits = new ArrayList<>();
		while( tokenizer.hasMoreTokens() ){
			strLimits.add( tokenizer.nextToken() );
		}

		return new KSTermini( strand, numFlexRes, strLimits );

	}


	public ArrayList<String> getWTSequence() {

		Molecule m = PDBFileReader.readPDBFile(params.getValue("PDBNAME"), null);
		ArrayList<String> flexibleRes = getFlexRes();
		int numPos = flexibleRes.size();
		ArrayList<String> wt = new ArrayList<>(); for( int pos = 0; pos < numPos; ++pos ) wt.add(null);

		for(int pos=0; pos<numPos; pos++) {
			Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
			String wtName = res.template.name;
			wt.set(pos, wtName);
		}

		wt.trimToSize();
		return wt;
	}


	public SearchProblem getSearchProblem( int strand, KSAllowedSeqs strandSeqs ) {

		String tmp = getParams().getValue("kStarOutputDir");
		if(tmp.equalsIgnoreCase("runName")) tmp = getParams().getValue("RUNNAME");

		String ematDir = tmp + File.separator + getParams().getValue("kStarEmatDir");
		ObjectIO.makeDir(ematDir, getParams().getBool("kStarDeleteEmatDir", false));
		String name = ematDir + File.separator + getParams().getValue("RUNNAME");

		String suffix = KSTermini.getTerminiString(strand);

		ArrayList<String[]> moveableStrands = strandSeqs.getMoveableStrandTermini();
		ArrayList<String[]> freeBBZones = strandSeqs.getFreeBBZoneTermini();
		DEEPerSettings dset = strandSeqs.getDEEPerSettings();

		System.out.println("CREATING SEARCH PROBLEM.  NAME: "+name);

		return new SearchProblem( name+"."+suffix, getParams().getValue("PDBNAME"), 
				strandSeqs.getFlexRes(), 
				strandSeqs.getAllowedAAs(),
				getParams().getBool("AddWT"), 
				getParams().getBool("AddWT", true), 
				getParams().getBool("doMinimize", false),
				new EPICSettings(params),
				getParams().getBool("UseTupExp", false),
				new LUTESettings(params),
				dset, moveableStrands, freeBBZones,
				getParams().getBool("useEllipses", false),
				getParams().getBool("useERef", false),
				getParams().getBool("AddResEntropy", false),
				getParams().getBool("addWTRots", false),
				getStrandLimits(strand),
				getParams().getBool("useVoxelG", false),
                                new ArrayList<>()
				);
	}


	ArrayList<String> getFlexResByStrand( int strand ) {

		if( strand != KSTermini.COMPLEX && strand != KSTermini.PROTEIN && strand != KSTermini.LIGAND )
			throw new RuntimeException("ERROR: specified strand " + strand + " is invalid");

		ArrayList<String> flexResList = new ArrayList<>();

		String resListString = params.getValue("strandMut"+strand);
		StringTokenizer tokenizer = new StringTokenizer(resListString);

		while(tokenizer.hasMoreTokens()){
			flexResList.add( tokenizer.nextToken() );
		}

		return flexResList;
	}


	public void verifyStrandsMutuallyExclusive() {
		// make sure that strands are mutually exclusive. assuming only two strands for now...
		KSTermini s0 = getStrandLimits(0);
		KSTermini s1 = getStrandLimits(1);

		// there is a C common to both integer sets
		// x1 <= C <= x2
		// y1 <= C <= y2
		// x1 <= y2 && y1 <= x2
		if(s0.getTerminusBegin() <= s1.getTerminusEnd() && s1.getTerminusBegin() <= s0.getTerminusEnd())
			throw new RuntimeException("ERROR: strand0 overlaps with strand1. Please fix strand termini.");
	}


	public KSAllowedSeqs getAllowedSequences(int strand, KSAllowedSeqs complexSeqs) {

		KSTermini limits = getStrandLimits(strand);

		if( complexSeqs == null ) {
			// this is only true when we have not calculated the sequences for the complex

			ArrayList<String> flexRes = getFlexRes();
			ArrayList<ArrayList<String>> allowedAAs = getAllowedAAs();

			if(flexRes.size() != allowedAAs.size()){
				throw new RuntimeException("ERROR: Number of flexible positions different in flexible residue "
						+ "and allowed AA type parameters!");
			}

			int numMutations = params.getInt("NUMMUTATIONS", 1);
                        boolean allowLessMut = params.getBool("ALLOWLESSMUTATIONS",false);

			complexSeqs = new KSAllowedSeqs(strand, limits, setupDEEPer(), 
					freeBBZoneTermini(limits), moveableStrandTermini(limits), flexRes, 
					allowedAAs, getWTSequence(), getParams().getBool("addWT"), 
                                        numMutations, allowLessMut);

			if( !complexSeqs.containsWTSeq() ) {
				System.out.println("WARNING: allowed sequences does not contain the wild-type sequence: " + 
						KSAbstract.list1D2String(complexSeqs.getWTSeq(), " ") + "\n");
			}

			// if this condition is true, then only the wild type sequence is returned
			if(numMutations > 0 && complexSeqs.getNumSeqs() == 1 && complexSeqs.containsWTSeq())
				throw new RuntimeException("ERROR: cannot generate any sequences "
						+ "for NUMMUTATIONS=" + numMutations + " mutation(s). "
						+ "Change the value of NUMMUTATIONS parameter.");
		}

		if(strand == KSTermini.COMPLEX)
			return complexSeqs;

		// get index values of strand limits in flexRes
		@SuppressWarnings("unchecked")
		ArrayList<String> flexRes = (ArrayList<String>)ObjectIO.deepCopy(complexSeqs.getFlexRes());
		@SuppressWarnings("unchecked")
		ArrayList<ArrayList<String>> allowedAAs = (ArrayList<ArrayList<String>>)ObjectIO.deepCopy(complexSeqs.getAllowedAAs());
		int lb = -1, ub = -1;
		KSTermini compressedLimits = getCompressedStrandLimits(strand, complexSeqs.getFlexRes());
		lb = flexRes.indexOf( String.valueOf(compressedLimits.getTerminusBegin()) );
		ub = flexRes.indexOf( String.valueOf(compressedLimits.getTerminusEnd()) ) + 1;

		// filter allowedSeqs for protein and ligand strands
		filterResiduesByStrand( strand, flexRes, allowedAAs );

		KSAllowedSeqs strandSeqs = new KSAllowedSeqs(strand, getStrandLimits(strand), setupDEEPer(strand), 
				freeBBZoneTermini(limits), moveableStrandTermini(limits), flexRes, complexSeqs, allowedAAs, lb, ub);

		return strandSeqs;
	}


	KSTermini getCompressedStrandLimits( int strand, ArrayList<String> flexRes ) {

		KSTermini strandLimits = getStrandLimits(strand);
		ArrayList<String> strandResNums = getFlexResByStrand(strand);

		@SuppressWarnings("unchecked")
		ArrayList<String> flexRes2 = (ArrayList<String>)ObjectIO.deepCopy(flexRes);

		for( int it = 0; it < flexRes2.size(); ) {
			if( !strandLimits.contains( Integer.parseInt(flexRes2.get(it)) ) ) {

				String removed = flexRes2.remove(it);

				if( strandResNums.contains( removed ) )
					throw new RuntimeException("ERROR: in strand" + strand + 
							", flexible residue " + removed + " exceeds strand limits");
			}
			else ++it;
		}

		ArrayList<String> limits = new ArrayList<>();
		limits.add(flexRes2.get(0));
		limits.add(flexRes2.get(flexRes2.size()-1));
		KSTermini ans = new KSTermini(strand, flexRes2.size(), limits);

		return ans;
	}


	void filterResiduesByStrand( int strand, ArrayList<String> flexRes, 
			ArrayList<ArrayList<String>> allowedAAs ) {

		KSTermini strandLimits = getStrandLimits(strand);
		ArrayList<String> strandResNums = getFlexResByStrand(strand);

		for( int it = 0; it < flexRes.size(); ) {
			if( !strandLimits.contains( Integer.parseInt(flexRes.get(it)) ) ) {

				String removed = flexRes.remove(it);

				if( strandResNums.contains( removed ) )
					throw new RuntimeException("ERROR: in strand" + strand + 
							", flexible residue " + removed + " exceeds strand limits");

				allowedAAs.remove(it);
			}
			else ++it;
		}

		if(flexRes.size() != allowedAAs.size()){
			throw new RuntimeException("ERROR: Number of flexible positions different in flexible residue "
					+ "and allowed AA type parameters!");
		}
	}
}
