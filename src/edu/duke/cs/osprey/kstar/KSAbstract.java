package edu.duke.cs.osprey.kstar;

import java.io.File;
import java.math.BigInteger;
import java.net.InetAddress;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

import edu.duke.cs.osprey.confspace.AllowedSeqs;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.confspace.SearchProblem.MatrixType;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;


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
	protected HashMap<Integer, HashMap<ArrayList<String>, PFAbstract>> strand2PFs = new HashMap<>();

	private ArrayList<SearchProblem> allSPs = new ArrayList<>();

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

		strand2PFs.put(Strand.COMPLEX, new HashMap<ArrayList<String>, PFAbstract>());
		strand2PFs.put(Strand.PROTEIN, new HashMap<ArrayList<String>, PFAbstract>());
		strand2PFs.put(Strand.LIGAND, new HashMap<ArrayList<String>, PFAbstract>());
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

	public String getSearchProblemName(int strand, ArrayList<String> seq) {
		return getEMATdir() + 
				File.separator + 
				getRunName() +"."+ 
				Strand.getStrandString(strand) + "." + 
				arrayList1D2String(seq, ".");
	}

	public String getEmatName(int strand, ArrayList<String> seq, SearchProblem.MatrixType type) {
		return getSearchProblemName(strand, seq) + "." + type.name() + ".dat";
	}

	protected void printSequences() {
		System.out.println("\nPreparing to compute K* for the following sequences:");
		int i = 0;
		for(ArrayList<String> al : strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqList()) {
			System.out.println(i++ + "\t" + arrayList1D2String(al, " "));
		}
		System.out.println();
	}

	private synchronized void addSPToGlobalList(int strand, SearchProblem sp) {
		int arrayPos = strand == Strand.COMPLEX ? 0 : Math.max(allSPs.size()-1, 0);
		allSPs.add(arrayPos, sp);

		if(allSPs.size() >= ThreadParallelism.getNumThreads()*2) {
			// create energy matrices and clear list to avoid running out of memory
			loadEnergyMatrices();
		}
	}

	private synchronized void addSPToLocalMap(int strand, SearchProblem sp, HashMap<Integer, SearchProblem> map) {
		map.put(strand, sp);
	}

	private synchronized void addPFToGlobalMap(int strand, ArrayList<String> seq, PFAbstract pf) {
		strand2PFs.get(strand).put(seq, pf);
	}

	private synchronized void addPFToLocalMap(int strand, PFAbstract pf, HashMap<Integer, PFAbstract> map) {
		map.put(strand, pf);
	}

	private void loadEnergyMatrices() {

		try {

			ForkJoinPool forkJoinPool = new ForkJoinPool(ThreadParallelism.getNumThreads());
			forkJoinPool.submit(() -> allSPs.parallelStream().forEach(sp -> {
				sp.loadEnergyMatrix();
			})).get();

			allSPs.clear(); 
		} 

		catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		} 
	}

	protected HashMap<Integer, PFAbstract> createPartitionFunctionsForSeq(int i) {

		HashMap<Integer, PFAbstract> ans = new HashMap<Integer, PFAbstract>();

		ArrayList<Integer> strands = new ArrayList<>();
		strands.add(Strand.COMPLEX);
		strands.add(Strand.PROTEIN);
		strands.add(Strand.LIGAND);

		strands.parallelStream().forEach((strand) -> {

			AllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);

			ArrayList<String> seq = strandSeqs.getStrandSeq(i);
			if( strand2PFs.get(strand).get(seq) == null ) {

				ArrayList<ArrayList<String>> allowedAAs = arrayList1D2ListOfLists(seq);
				ArrayList<Integer> pos = new ArrayList<>(); for( int j = 0; j < seq.size(); ++j ) pos.add(j);
				ArrayList<String> flexibleRes = strand2AllSearchProblem.get(strand).getFlexibleResiduePositions(seq, pos);
				ArrayList<String[]> moveableStrands = strand2AllowedSeqs.get(strand).getMoveableStrandTermini();
				ArrayList<String[]> freeBBZones = strand2AllowedSeqs.get(strand).getFreeBBZoneTermini();
				DEEPerSettings dset = strand2AllowedSeqs.get(strand).getDEEPerSettings();

				// create searchproblem
				SearchProblem seqSearchProblem = new SearchProblem( 
						getSearchProblemName(strand, seq), 
						pdbName, 
						flexibleRes, 
						allowedAAs, 
						true, 
						true,
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
						PFAbstract.getImplementation(),
						seq, 
						cfp, 
						seqSearchProblem, 
						null, 
						dset, 
						moveableStrands, 
						freeBBZones,
						EW+I0);

				// put partition function in global map
				addPFToGlobalMap(strand, seq, pf);

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
		for(int strand : strands) {

			if(!ans.keySet().contains(strand)) {

				AllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);
				ArrayList<String> seq = strandSeqs.getStrandSeq(i);

				PFAbstract pf = strand2PFs.get(strand).get(seq);

				ans.put(strand, pf);
			}
		}

		if(ans.size() != 3)
			throw new RuntimeException("ERROR: returned map must contain three different partition functions");

		return ans;
	}

	protected HashMap<Integer, SearchProblem> createSearchProblemsForSeq(int i) {
		// used to precompute energy matrices
		HashMap<Integer, SearchProblem> ans = new HashMap<Integer, SearchProblem>();

		ArrayList<Integer> strands = new ArrayList<>();
		strands.add(Strand.COMPLEX);
		strands.add(Strand.PROTEIN);
		strands.add(Strand.LIGAND);

		strands.parallelStream().forEach((strand) -> {

			AllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);

			ArrayList<String> seq = strandSeqs.getStrandSeq(i);
			if( !new File(getEmatName(strand, seq, MatrixType.EMAT)).exists() ) {

				ArrayList<ArrayList<String>> allowedAAs = arrayList1D2ListOfLists(seq);
				ArrayList<Integer> pos = new ArrayList<>(); for( int j = 0; j < seq.size(); ++j ) pos.add(j);
				ArrayList<String> flexibleRes = strand2AllSearchProblem.get(strand).getFlexibleResiduePositions(seq, pos);
				ArrayList<String[]> moveableStrands = strand2AllowedSeqs.get(strand).getMoveableStrandTermini();
				ArrayList<String[]> freeBBZones = strand2AllowedSeqs.get(strand).getFreeBBZoneTermini();
				DEEPerSettings dset = strand2AllowedSeqs.get(strand).getDEEPerSettings();

				// create searchproblem
				SearchProblem seqSearchProblem = new SearchProblem( 
						getSearchProblemName(strand, seq), 
						pdbName, 
						flexibleRes, 
						allowedAAs, 
						true, 
						true,
						useEPIC,
						new EPICSettings(cfp.getParams()),
						useTupExp,
						dset, 
						moveableStrands, 
						freeBBZones,
						useEllipses,
						useERef,
						addResEntropy);

				// synchronized
				addSPToLocalMap(strand, seqSearchProblem, ans);
			}		
		});

		return ans;
	}


	protected void createEnergyMatrices() {

		System.out.println("\nCreating all energy matrices\n");

		try {

			ForkJoinPool forkJoinPool = new ForkJoinPool(ThreadParallelism.getNumThreads());
			forkJoinPool.submit(() ->

			IntStream.range(0, strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs()).parallel().forEach(i -> {

				System.out.println("\nCreating search problem for sequence " + 
						i + "/" + strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs() + "\n");

				HashMap<Integer, SearchProblem> map = createSearchProblemsForSeq(i);

				// put partition function in list, so we can parallelize energy matrix computation
				for(int strand : map.keySet()) {
					addSPToGlobalList(strand, map.get(strand));
				}

			})).get();

			// create last of the energy matrices
			loadEnergyMatrices();

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
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

		for(HashMap<ArrayList<String>, PFAbstract> map : strand2PFs.values()) {
			for(PFAbstract pf : map.values()) {
				ans = ans.add(pf.getNumMinimizedConfs());
				ans = ans.add(pf.getPhaseOneMinimizedConfs());
			}
		}

		return ans;
	}

	protected BigInteger countTotNumConfs() {
		BigInteger ans = BigInteger.ZERO;

		for(HashMap<ArrayList<String>, PFAbstract> map : strand2PFs.values()) {
			for(PFAbstract pf : map.values()) {
				ans = ans.add(pf.getNumInitialUnPrunedConfs());
			}
		}

		return ans;
	}
}
