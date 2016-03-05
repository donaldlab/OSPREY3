package edu.duke.cs.osprey.kstar;

import java.io.File;
import java.math.BigInteger;
import java.net.InetAddress;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ForkJoinPool;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.kstar.pfunction.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunction.PFFactory;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import edu.duke.cs.osprey.tools.ObjectIO;


/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public abstract class KSAbstract implements KSInterface {

	protected ConfigFileParser cfp = null;
	protected ArrayList<String> wtSeq = null;
	protected KSCalc wtKSCalc = null;
	protected String outFileName = null;
	protected String ematDir = null;
	protected String runName = null;

	protected HashMap<Integer, SearchProblem> strand2AllSearchProblem = new HashMap<Integer, SearchProblem>();
	protected HashMap<Integer, AllowedSeqs> strand2AllowedSeqs = null;
	protected HashMap<Boolean, HashMap<Integer, HashMap<ArrayList<String>, PFAbstract>>> contSCFlex2PFs = new HashMap<>();

	private ArrayList<SearchProblem> tempSPs = new ArrayList<>();
	protected HashSet<String> allSPNames = new HashSet<>();

	double EW;
	double I0;
	private String pdbName;
	private boolean useEPIC;
	private boolean useTupExp;
	private boolean useEllipses;
	private boolean useERef;
	private boolean addResEntropy;


	public KSAbstract( ConfigFileParser cfp ) {

		this.cfp = cfp;

		EW = cfp.getParams().getDouble("Ew",0);
		I0 = cfp.getParams().getBool("imindee",false) ? cfp.getParams().getDouble("Ival",5) : 0;
		pdbName = cfp.getParams().getValue("PDBNAME");
		useEPIC = cfp.getParams().getBool("UseEPIC");
		useTupExp = cfp.getParams().getBool("UseTupExp");
		useEllipses = cfp.getParams().getBool("useEllipses");
		useERef = cfp.getParams().getBool("useERef");
		addResEntropy = cfp.getParams().getBool("AddResEntropy");

		// populate map
		boolean[] contSCFlexVals = { true, false };
		int[] strands = { Strand.COMPLEX, Strand.PROTEIN, Strand.LIGAND };

		for( boolean contSCFlex : contSCFlexVals ) {

			contSCFlex2PFs.put(contSCFlex, new HashMap<Integer, HashMap<ArrayList<String>, PFAbstract>>());

			for( int strand : strands )
				contSCFlex2PFs.get(contSCFlex).put(strand, new HashMap<ArrayList<String>, PFAbstract>());
		}
	}

	protected void createEmatDir() {
		if( cfp.getParams().getBool("deleteematdir", false) )
			ObjectIO.deleteDir(getEMATdir());

		if( !new File(getEMATdir()).exists() )
			ObjectIO.makeDir(getEMATdir(), false);
	}

	public static String arrayList1D2String(ArrayList<String> seq, String separator) {
		StringBuilder ans = new StringBuilder();
		for( int i = 0; i < seq.size(); ++i ) {
			if(i > 0) ans.append(separator);
			ans.append(seq.get(i));
		}

		return ans.toString();
	}

	public static ArrayList<ArrayList<String>> arrayList1D2ListOfLists(ArrayList<String> list) {
		ArrayList<ArrayList<String>> ans = new ArrayList<>();
		for( String s : list ) {
			ArrayList<String> newList = new ArrayList<>();
			newList.add(s);
			ans.add(newList);
		}

		return ans;
	}

	public String getSearchProblemName(boolean contSCFlex, int strand) {

		String flexibility = contSCFlex == true ? "min" : "rig";

		return getEMATdir() + 
				File.separator + 
				getRunName() + "." + 
				flexibility + "." + 
				Strand.getStrandString(strand);
	}

	public String getSearchProblemName(boolean contSCFlex, int strand, ArrayList<String> seq) {
		return getSearchProblemName(contSCFlex, strand) + "." + arrayList1D2String(seq, ".");
	}

	public String getEmatName(boolean contSCFlex, int strand, ArrayList<String> seq, SearchProblem.MatrixType type) {
		return getSearchProblemName(contSCFlex, strand, seq) + "." + type.name() + ".dat";
	}

	protected void printSequences() {
		System.out.println("\nPreparing to compute K* for the following sequences:");
		int i = 0;
		for(ArrayList<String> al : strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqList()) {
			System.out.println(i++ + "\t" + arrayList1D2String(al, " "));
		}
		System.out.println();
	}

	protected synchronized void addSPToTmpList(int strand, SearchProblem sp) {
		int arrayPos = strand == Strand.COMPLEX ? 0 : Math.max(tempSPs.size()-1, 0);
		tempSPs.add(arrayPos, sp);

		if(tempSPs.size() >= ThreadParallelism.getNumThreads()*2) {
			// create energy matrices and clear list to avoid running out of memory
			loadEnergyMatrices();
		}
	}

	protected synchronized void addSPToLocalMap(int strand, SearchProblem sp, HashMap<Integer, SearchProblem> map) {
		map.put(strand, sp);
	}

	private synchronized void addPFToGlobalMap(boolean contSCFlex, int strand, ArrayList<String> seq, PFAbstract pf) {
		contSCFlex2PFs.get(contSCFlex).get(strand).put(seq, pf);
	}

	private synchronized void addPFToLocalMap(int strand, PFAbstract pf, HashMap<Integer, PFAbstract> map) {
		map.put(strand, pf);
	}

	protected void loadEnergyMatrices() {

		try {

			ForkJoinPool forkJoinPool = new ForkJoinPool(ThreadParallelism.getNumThreads());
			forkJoinPool.submit(() -> tempSPs.parallelStream().forEach(sp -> {

				sp.loadEnergyMatrix();

			})).get();

			tempSPs.clear(); 
		} 

		catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		} 
	}

	protected HashMap<Integer, PFAbstract> createPartitionFunctionsForSeq(String[][] seqs, 
			boolean[] contSCFlexVals, String[] pfImplVals) {

		HashMap<Integer, PFAbstract> ans = new HashMap<Integer, PFAbstract>();

		int[] strands = { Strand.COMPLEX, Strand.PROTEIN, Strand.LIGAND };
		ArrayList<Integer> indexes = new ArrayList<>();
		for(int i = 0; i < strands.length; i++) indexes.add(i);

		indexes.parallelStream().forEach((index) -> {

			int strand = strands[index];
			boolean contSCFlex = contSCFlexVals[index];
			String pfImpl = pfImplVals[index];

			ArrayList<String> seq = new ArrayList<String>(Arrays.asList(seqs[strand]));
			
			if( contSCFlex2PFs.get(contSCFlex).get(strand).get(seq) == null ) {

				ArrayList<ArrayList<String>> allowedAAs = arrayList1D2ListOfLists(seq);
				// ArrayList<Integer> pos = new ArrayList<>(); for( int j = 0; j < seq.size(); ++j ) pos.add(j);
				// ArrayList<String> flexibleRes = strand2AllSearchProblem.get(strand).getFlexibleResiduePositions(seq, pos);
				ArrayList<String> flexibleRes = strand2AllowedSeqs.get(strand).getFlexRes(seq);
				ArrayList<String[]> moveableStrands = strand2AllowedSeqs.get(strand).getMoveableStrandTermini();
				ArrayList<String[]> freeBBZones = strand2AllowedSeqs.get(strand).getFreeBBZoneTermini();
				DEEPerSettings dset = strand2AllowedSeqs.get(strand).getDEEPerSettings();

				// create searchproblem
				SearchProblem seqSearchProblem = new SearchProblem( 
						getSearchProblemName(contSCFlex, strand, seq), 
						pdbName, 
						flexibleRes, 
						allowedAAs, 
						false, 
						contSCFlex,
						useEPIC,
						new EPICSettings(cfp.getParams()),
						useTupExp,
						dset, 
						moveableStrands, 
						freeBBZones,
						useEllipses,
						useERef,
						addResEntropy);

				// create partition function
				PFAbstract pf = PFFactory.getPartitionFunction( 
						pfImpl,
						seq, 
						cfp, 
						seqSearchProblem, 
						null, 
						dset, 
						moveableStrands, 
						freeBBZones,
						EW+I0);

				// put partition function in global map
				addPFToGlobalMap(contSCFlex, strand, seq, pf);

				// put in local map
				addPFToLocalMap(strand, pf, ans);

				// get energy matrix
				pf.getSearchProblem().loadEnergyMatrix();

				// get pruning matrix
				pf.getPruningControl().prune();

				// initialize conf counts for K*
				pf.setNumUnPrunedConfs();
				pf.setNumPrunedConfsByDEE();
			}		
		});

		// get pfs that were already in global map
		for( int index : indexes ) {

			int strand = strands[index];
			boolean contSCFlex = contSCFlexVals[index];

			if(!ans.keySet().contains(strand)) {

				ArrayList<String> seq = new ArrayList<String>(Arrays.asList(seqs[strand]));
				
				PFAbstract pf = contSCFlex2PFs.get(contSCFlex).get(strand).get(seq);

				ans.put(strand, pf);
			}
		}

		if(ans.size() != 3)
			throw new RuntimeException("ERROR: returned map must contain three different partition functions");

		return ans;
	}


	protected boolean createSP( String name ) {

		synchronized( allSPNames ) {

			if(!allSPNames.contains(name)) {

				allSPNames.add(name);

				return true;
			}

			return false;
		}
	}


	protected SearchProblem createSingleSequenceSearchProblem( boolean contSCFlex, int strand, ArrayList<String> seq ) {

		ArrayList<ArrayList<String>> allowedAAs = arrayList1D2ListOfLists(seq);
		// ArrayList<Integer> pos = new ArrayList<>(); for( int j = 0; j < seq.size(); ++j ) pos.add(j);
		// ArrayList<String> flexibleRes = strand2AllSearchProblem.get(strand).getFlexibleResiduePositions(seq, pos);
		ArrayList<String> flexibleRes = strand2AllowedSeqs.get(strand).getFlexRes(seq);
		ArrayList<String[]> moveableStrands = strand2AllowedSeqs.get(strand).getMoveableStrandTermini();
		ArrayList<String[]> freeBBZones = strand2AllowedSeqs.get(strand).getFreeBBZoneTermini();
		DEEPerSettings dset = strand2AllowedSeqs.get(strand).getDEEPerSettings();

		// create searchproblem
		SearchProblem seqSearchProblem = new SearchProblem( 
				getSearchProblemName(contSCFlex, strand, seq), 
				pdbName, 
				flexibleRes, 
				allowedAAs, 
				false, 
				contSCFlex,
				useEPIC,
				new EPICSettings(cfp.getParams()),
				useTupExp,
				dset, 
				moveableStrands, 
				freeBBZones,
				useEllipses,
				useERef,
				addResEntropy);

		return seqSearchProblem;
	}


	protected SearchProblem createAllSequenceSearchProblem( boolean contSCFlex, int strand ) {

		ArrayList<ArrayList<String>> allowedAAs = strand2AllowedSeqs.get(strand).getAllowedAAs();
		ArrayList<String> flexibleRes = strand2AllowedSeqs.get(strand).getFlexRes();
		ArrayList<String[]> moveableStrands = strand2AllowedSeqs.get(strand).getMoveableStrandTermini();
		ArrayList<String[]> freeBBZones = strand2AllowedSeqs.get(strand).getFreeBBZoneTermini();
		DEEPerSettings dset = strand2AllowedSeqs.get(strand).getDEEPerSettings();

		// create searchproblem
		SearchProblem allSeqSearchProblem = new SearchProblem( 
				getSearchProblemName(contSCFlex, strand), 
				pdbName, 
				flexibleRes, 
				allowedAAs, 
				true, 
				contSCFlex,
				useEPIC,
				new EPICSettings(cfp.getParams()),
				useTupExp,
				dset, 
				moveableStrands, 
				freeBBZones,
				useEllipses,
				useERef,
				addResEntropy);

		return allSeqSearchProblem;
	}


	protected String getOputputFileName() {

		if(outFileName == null) {
			try {

				outFileName = getRunName() 
						+ "." + InetAddress.getLocalHost().getHostName()
						+ "." + getKSMethod()
						+ "." + cfp.getParams().getValue("pFuncMethod", "1npcpmcache")
						+ ".txt";

			} catch (Exception e) {
				System.out.println(e.getMessage());
				e.printStackTrace();
				System.exit(1);
			}
		}

		return outFileName;
	}

	protected boolean getConcurrentNodeExpansion() {

		String expansionType = cfp.getParams().getValue("KStarExpansion", "serial");

		return expansionType.equalsIgnoreCase("serial") ? false : true;
	}

	protected String getRunName() {
		if(runName == null) runName = cfp.getParams().getValue("runName");
		return runName;
	}

	protected String getEMATdir() {
		if(ematDir == null) ematDir = cfp.getParams().getValue("ematdir", "emat");
		return ematDir;
	}

	protected ArrayList<String> getWTSeq() {
		if( wtSeq == null ) wtSeq = cfp.getWTSequence();
		return wtSeq;
	}

	protected KSCalc getWTKSCalc() {
		return wtKSCalc;
	}

	protected boolean isWT( KSCalc mutSeq ) {
		return mutSeq.getPF(Strand.COMPLEX).getSequence().equals(getWTSeq());
	}

	protected BigInteger countMinimizedConfs() {

		BigInteger ans = BigInteger.ZERO;

		for(HashMap<Integer, HashMap<ArrayList<String>, PFAbstract>> strand : contSCFlex2PFs.values()) {

			for(HashMap<ArrayList<String>, PFAbstract> seq : strand.values()) {

				for(PFAbstract pf : seq.values()) {

					ans = ans.add(pf.getNumMinimizedConfs());
				}
			}
		}

		return ans;
	}

	protected BigInteger countTotNumConfs() {

		BigInteger ans = BigInteger.ZERO;

		for(HashMap<Integer, HashMap<ArrayList<String>, PFAbstract>> strand : contSCFlex2PFs.values()) {

			for(HashMap<ArrayList<String>, PFAbstract> seq : strand.values()) {

				for(PFAbstract pf : seq.values()) {
					ans = ans.add(pf.getNumInitialUnPrunedConfs());
				}
			}
		}

		return ans;
	}
	
	protected String[][] getStrandStringsAtPos(int i) {
		
		String[][] ans = new String[3][];
		
		ans[Strand.COMPLEX] = (String[]) strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqAtPos(i).toArray(new String[0]);
		ans[Strand.PROTEIN] = (String[]) strand2AllowedSeqs.get(Strand.PROTEIN).getStrandSeqAtPos(i).toArray(new String[0]);
		ans[Strand.LIGAND] = (String[]) strand2AllowedSeqs.get(Strand.LIGAND).getStrandSeqAtPos(i).toArray(new String[0]);
		
		return ans;
	}
}
