package edu.duke.cs.osprey.kstar;

import java.io.File;
import java.math.BigInteger;
import java.net.InetAddress;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;
import java.util.concurrent.ConcurrentHashMap;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.RunState;
import edu.duke.cs.osprey.kstar.pfunc.PFFactory;
import edu.duke.cs.osprey.pruning.PruningControl;
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
	protected String outFilePath = null;
	protected String ematDir = null;
	protected String checkPointDir = null;
	protected String runName = null;
	protected String checkPointFilePath = null;

	protected HashMap<Integer, AllowedSeqs> strand2AllowedSeqs = new HashMap<>();
	protected ConcurrentHashMap<String, SearchProblem> name2SP = new ConcurrentHashMap<>();
	protected ConcurrentHashMap<String, PFAbstract> name2PF = new ConcurrentHashMap<>();

	protected double EW;
	protected double I0;
	private String pdbName;
	private boolean useEPIC;
	private boolean useTupExp;
	private boolean useEllipses;
	private boolean useERef;
	private boolean addResEntropy;

	private static double pRatioLBT = 0.25;
	private static double pRatioUBT = 0.95;
	protected boolean prunedSingleSeqs = false;
	public static boolean preLoadPFs = false;
	public static boolean refinePruning = false;
	public static boolean doCheckPoint = false;
	protected static long checkpointInterval = 50000;

	public static long begin = 0; // beginning time

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
	}


	protected abstract void preLoadPFs(ArrayList<Boolean> contSCFlexVals);


	public static void setCheckPointInterval( long interval ) {
		if(interval <= 0) interval = 1;
		checkpointInterval = interval;
	}


	public void createEmats(ArrayList<Boolean> contSCFlexVals) {
		// for now, only the pan seqSPs have energy matrices in the file system
		preparePanSeqSPs(contSCFlexVals);

		// for benchmarking, prepare these ahead of time
		// only single sequences are loaded and pruned this time
		if(preLoadPFs) preLoadPFs(contSCFlexVals);
	}


	protected void createCheckPointDir() {
		if( !new File(getCheckPointDir()).exists() )
			ObjectIO.makeDir(getCheckPointDir(), false);
	}


	protected void createEmatDir() {
		ObjectIO.makeDir(getEMATdir(), cfp.getParams().getBool("deleteematdir", false));
	}


	public static String list1D2String(ArrayList<String> seq, String separator) {
		StringBuilder ans = new StringBuilder();
		for( int i = 0; i < seq.size(); ++i ) {
			if(i > 0) ans.append(separator);
			ans.append(seq.get(i));
		}

		return ans.toString();
	}


	public static ArrayList<String> file2List( String path ) {
		ArrayList<String> ans = new ArrayList<String>();

		try {
			if( !new File(path).exists() ) return ans;

			Scanner s = new Scanner(new File(path));
			while (s.hasNextLine()) ans.add(s.nextLine());
			s.close();

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

		return ans;
	}


	public static ArrayList<ArrayList<String>> list1D2ListOfLists(ArrayList<String> list) {
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


	private String getCheckPointName(boolean contSCFlex, int strand) {

		String flexibility = contSCFlex == true ? "min" : "rig";

		return getCheckPointDir() + 
				File.separator + 
				getRunName() + "." + 
				flexibility + "." + 
				Strand.getStrandString(strand);
	}


	public String getSearchProblemName(boolean contSCFlex, int strand, String pfImpl, ArrayList<String> seq) {
		return getSearchProblemName(contSCFlex, strand) + "." + pfImpl + "." + list1D2String(seq, ".");
	}


	public String getCheckPointName(boolean contSCFlex, int strand, String pfImpl, ArrayList<String> seq) {
		return getCheckPointName(contSCFlex, strand) + "." + pfImpl + "." + list1D2String(seq, ".") +".checkpoint";
	}


	protected void printSequences() {
		System.out.println("\nPreparing to compute K* for the following sequences:");
		int i = 0;
		for(ArrayList<String> al : strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqList()) {
			System.out.println(i++ + "\t" + list1D2String(al, " "));
		}
		System.out.println();
	}


	protected void loadAndPruneMatricesFromSPMap() {

		try {

			System.out.println("\nCreating and pruning pan energy matrices\n");
			long begin = System.currentTimeMillis();

			name2SP.keySet().parallelStream().forEach(key -> {

				SearchProblem sp = name2SP.get(key);

				// single seq matrices created using the fast construction 
				// method are already pruned according to the pruning window
				if(sp.getEnergyMatrix() == null) {
					PruningControl pc = cfp.getPruningControl(sp, EW+I0, useEPIC, useTupExp); 
					sp.loadEnergyMatrix(sp.getMatrixType());
					pc.prune();
				}

				// for K* we will guard against under-pruning, since this increases
				// the length of our calculation
				if(sp.isSingleSeq() && KSAbstract.refinePruning) {
					refinePruningInterval(sp);
				}
			});

			System.out.println("\nFinished creating and pruning energy matrices");
			System.out.println("Running time: " + (System.currentTimeMillis()-begin)/1000 + " seconds\n");
		} 

		catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		} 
	}


	protected void loadAndPruneMatricesFromPFMap() {

		try {
			System.out.println("\nCreating and pruning pan energy matrices\n");
			long begin = System.currentTimeMillis();

			name2PF.keySet().parallelStream().forEach(key -> {

				SearchProblem sp = name2PF.get(key).getSearchProblem();

				// single seq matrices created using the fast construction 
				// method are already pruned according to the pruning window
				if(sp.getEnergyMatrix() == null) {
					PruningControl pc = cfp.getPruningControl(sp, EW+I0, useEPIC, useTupExp); 
					sp.loadEnergyMatrix(sp.getMatrixType());
					pc.prune();
				}

				// for K* we will guard against under-pruning, since this increases
				// the length of our calculation
				if(sp.isSingleSeq() && KSAbstract.refinePruning) {
					refinePruningInterval(sp);
				}
			});

			System.out.println("\nFinished creating and pruning energy matrices");
			System.out.println("Running time: " + (System.currentTimeMillis()-begin)/1000 + " seconds\n");
		} 

		catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		} 
	}


	protected void refinePruningInterval( SearchProblem sp ) {

		if(sp.numConfs(false).compareTo(BigInteger.ZERO) == 0)
			return;

		// if current interval is good enough, return
		double r = sp.numConfs(true).doubleValue() / sp.numConfs(false).doubleValue();
		if(r >= pRatioLBT && r <= pRatioUBT) return;

		// most amount of pruning; ratio upper bound
		double l = 0.01;
		PruningControl pc = cfp.getPruningControl(sp, l, useEPIC, useTupExp); pc.prune();
		double lr = sp.numConfs(true).doubleValue() / sp.numConfs(false).doubleValue();

		// least amount of pruning; ratio lower bound
		double u = 100;
		pc = cfp.getPruningControl(sp, u, useEPIC, useTupExp); pc.prune();
		double ur = sp.numConfs(true).doubleValue() / sp.numConfs(false).doubleValue();

		double m = -1, mr = -1;

		// ratio cannot get smaller or bigger, respectively
		if(ur > pRatioUBT || lr < pRatioLBT) {
			pc = cfp.getPruningControl(sp, EW+I0, useEPIC, useTupExp); pc.prune();
			sp.pruneMat.setPruningInterval(EW+I0);
			return;
		}

		while( Math.abs(l-u) > 0.1 && (lr > pRatioUBT || ur < pRatioLBT)) {

			m = (l+u)/2;
			pc = cfp.getPruningControl(sp, m, useEPIC, useTupExp); pc.prune();
			mr = sp.numConfs(true).doubleValue() / sp.numConfs(false).doubleValue();
			sp.pruneMat.setPruningInterval(m);

			if(mr < pRatioLBT) {
				ur = mr;
				u = m;
			}

			else if(mr > pRatioUBT) {
				lr = mr;
				l = m;
			}

			else if(mr >= pRatioLBT && mr <= pRatioUBT) {
				// tada!
				System.out.println(m + "\t" + mr);
				return;
			}
		}

		// we failed. restore original pruning interval
		pc = cfp.getPruningControl(sp, EW+I0, useEPIC, useTupExp); pc.prune();
		sp.pruneMat.setPruningInterval(EW+I0);
		return;
	}


	protected PFAbstract createPF4Seq(boolean contSCFlex, int strand, ArrayList<String> seq, String pfImpl) {
		PFAbstract ans = null;

		String panSeqSPName = getSearchProblemName(contSCFlex, strand);
		SearchProblem panSeqSP = name2SP.get(panSeqSPName);

		String seqSPName = getSearchProblemName(contSCFlex, strand, pfImpl, seq);
		String cpName = getCheckPointName(contSCFlex, strand, pfImpl, seq);

		if( (ans = name2PF.get(seqSPName)) != null ) return ans;
		if( doCheckPoint && (ans = deSerializePF(seqSPName, cpName, contSCFlex)) != null ) return ans;

		ArrayList<Integer> flexResIndexes = strand2AllowedSeqs.get(strand).getFlexResIndexesFromSeq(seq);

		// create partition function
		ans = PFFactory.getPartitionFunction(pfImpl, strand, seq, flexResIndexes, cpName, seqSPName, cfp, panSeqSP);

		return ans;
	}


	protected ConcurrentHashMap<Integer, PFAbstract> createPFs4Seqs(ArrayList<ArrayList<String>> seqs, 
			ArrayList<Boolean> contSCFlexVals, ArrayList<String> pfImplVals) {

		ConcurrentHashMap<Integer, PFAbstract> ans = new ConcurrentHashMap<>();

		ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(Strand.COMPLEX, Strand.PROTEIN, Strand.LIGAND));
		ArrayList<Integer> indexes = new ArrayList<>();
		for(int i = 0; i < strands.size(); i++) indexes.add(i);

		indexes.parallelStream().forEach((index) -> {
			//for(int index = 0; index < strands.size(); ++index) {

			int strand = strands.get(index);
			boolean contSCFlex = contSCFlexVals.get(strand);
			String pfImpl = pfImplVals.get(strand);
			ArrayList<String> seq = seqs.get(strand);

			PFAbstract pf = createPF4Seq(contSCFlex, strand, seq, pfImpl);

			// put partition function in global map
			name2PF.put(pf.getSearchProblemName(), pf);

			// put in local map
			ans.put(strand, pf);

			// only continue if we have not already started computed the PF
			if( pf.getRunState() == RunState.NOTSTARTED ) {

				// get energy matrix
				if(pf.getSearchProblem().getEnergyMatrix() == null) {
					pf.getSearchProblem().loadEnergyMatrix(pf.getSearchProblem().getMatrixType());

					// get pruning matrix
					pf.getPruningControl(EW+I0).prune();
				}

				if(pf.getSearchProblem().numConfs(false).compareTo(BigInteger.ZERO) == 0) {
					// no conformations in search space, so this cannot give a valid
					// partition function
					pf.setEpsilonStatus(EApproxReached.NOT_POSSIBLE);

					System.out.println("\nWARNING: there are no valid conformations for sequence " + 
							KSAbstract.list1D2String(pf.getSequence(), " "));
				}

				else {	
					if(!prunedSingleSeqs && KSAbstract.refinePruning) {
						refinePruningInterval(pf.getSearchProblem());
					}

					// initialize conf counts for K*
					pf.setNumUnPruned();
					pf.setNumPruned();
				}
			}
			//}
		});

		if(ans.size() != 3)
			throw new RuntimeException("ERROR: returned map must contain three different partition functions");

		return ans;
	}


	public void removeFromMap(String spName, boolean sp, boolean pf) {
		if(sp) name2SP.remove(spName);
		if(pf) name2PF.remove(spName);
	}


	protected SearchProblem createPanSeqSP( boolean contSCFlex, int strand ) {

		ArrayList<ArrayList<String>> allowedAAs = AllowedSeqs.removePosFromAllowedAAs(strand2AllowedSeqs.get(strand).getAllowedAAs());
		ArrayList<String> flexibleRes = strand2AllowedSeqs.get(strand).getFlexRes();
		ArrayList<String[]> moveableStrands = strand2AllowedSeqs.get(strand).getMoveableStrandTermini();
		ArrayList<String[]> freeBBZones = strand2AllowedSeqs.get(strand).getFreeBBZoneTermini();
		DEEPerSettings dset = strand2AllowedSeqs.get(strand).getDEEPerSettings();

		// create searchproblem
		SearchProblem panSeqSP = new SearchProblem( 
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

		return panSeqSP;
	}


	protected void preparePanSeqSPs( ArrayList<Boolean> contSCFlexVals ) {

		ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(Strand.LIGAND, 
				Strand.PROTEIN, Strand.COMPLEX));

		for( boolean contSCFlex : contSCFlexVals ) {

			strands.parallelStream().forEach(strand -> {

				String spName = getSearchProblemName(contSCFlex, strand);

				if( !name2SP.containsKey(spName) ) {
					SearchProblem sp = createPanSeqSP(contSCFlex, strand);
					name2SP.put(sp.name, sp);
				}
			});	
		}

		loadAndPruneMatricesFromSPMap();
	}


	protected String getOputputFilePath() {

		if(outFilePath == null)
			outFilePath = getRunName() + "." + getSummaryFileExtension();

		return outFilePath;
	}


	protected String getCheckPointFilePath() {

		if(checkPointFilePath == null)
			checkPointFilePath = getCheckPointDir() + File.separator + getOputputFilePath();

		return checkPointFilePath;
	}


	protected String getSummaryFileExtension() {

		String ans = null;
		try {

			ans = InetAddress.getLocalHost().getHostName()
					+ "." + getKSMethod()
					+ "." + PFAbstract.getCFGImpl()
					+ ".txt";

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

		return ans;
	}


	protected boolean getConcurrentNodeExpansion() {

		String expansionType = cfp.getParams().getValue("KStarExpansion", "true");

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


	protected synchronized String getCheckPointDir() {
		if(checkPointDir == null) checkPointDir = cfp.getParams().getValue("checkPointDir", "checkpoint");
		return checkPointDir;
	}


	protected ArrayList<String> getWTSeq() {
		if( wtSeq == null ) wtSeq = strand2AllowedSeqs.get(Strand.COMPLEX).getWTSeq();
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
		for(PFAbstract pf : name2PF.values()) ans = ans.add(pf.getNumMinimized4Output());
		return ans;
	}


	protected BigInteger countTotNumConfs() {

		BigInteger ans = BigInteger.ZERO;
		for(PFAbstract pf : name2PF.values()) ans = ans.add(pf.getNumUnPruned());
		return ans;
	}


	protected void abortPFs() {
		name2PF.keySet().parallelStream().forEach(key -> {
			PFAbstract pf = name2PF.get(key);
			pf.abort(true);
		});
	}


	protected ArrayList<ArrayList<String>> getStrandStringsAtPos(int i) {

		ArrayList<ArrayList<String>> ans = new ArrayList<ArrayList<String>>(Arrays.asList(null, null, null));

		ans.set(Strand.COMPLEX, strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqAtPos(i));
		ans.set(Strand.PROTEIN, strand2AllowedSeqs.get(Strand.PROTEIN).getStrandSeqAtPos(i));
		ans.set(Strand.LIGAND, strand2AllowedSeqs.get(Strand.LIGAND).getStrandSeqAtPos(i));

		ans.trimToSize();
		return ans;
	}


	protected KSCalc computeWTCalc() {
		
		if( !doCheckPoint || !new File(getOputputFilePath()).exists() ) KSCalc.printSummaryHeader(getOputputFilePath());
		if( doCheckPoint && !new File(getCheckPointFilePath()).exists() ) KSCalc.printSummaryHeader(getCheckPointFilePath());
		
		// compute wt sequence for reference
		ArrayList<ArrayList<String>> strandSeqs = getStrandStringsAtPos(0);		
		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, true, true));
		ArrayList<String> pfImplVals = new ArrayList<String>(Arrays.asList(PFAbstract.getCFGImpl(), 
				PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl()));

		ConcurrentHashMap<Integer, PFAbstract> pfs = createPFs4Seqs(strandSeqs, contSCFlexVals, pfImplVals);
		KSCalc calc = new KSCalc(0, pfs);
		
		PFAbstract pf = calc.getPF(Strand.COMPLEX);
		if(doCheckPoint && getSeqsFromFile(getOputputFilePath()).contains(pf.getSequence())) {
			// we have previously computed the sequence
			return calc;
		}
	
		calc.run(calc, false, true);

		// protein and ligand must reach epsilon, regardless of checkppoint
		for( int strand : Arrays.asList(Strand.LIGAND, Strand.PROTEIN, Strand.COMPLEX) ) {
			pf = calc.getPF(strand);

			if( (strand != Strand.COMPLEX && pf.getEpsilonStatus() != EApproxReached.TRUE) ||
					(strand == Strand.COMPLEX && pf.getEpsilonStatus() == EApproxReached.NOT_POSSIBLE) ) {
				
				throw new RuntimeException("ERROR: could not compute the wild-type sequence "
						+ KSAbstract.list1D2String(pf.getSequence(), " ") + " to an epsilon value of "
						+ PFAbstract.targetEpsilon + ". Resolve any clashes involving the flexible residues, "
						+"or increase the value of epsilon." );
			}
		}
		
		if( doCheckPoint ) {
			pf = calc.getPF(Strand.COMPLEX);
			name2PF.remove(pf.getSearchProblemName());
			calc.deleteSeqFromFile( pf.getSequence(), getCheckPointFilePath() );
			calc.serializePFs();
		}
		
		if( calc.getEpsilonStatus() == EApproxReached.TRUE ) {
			calc.deleteCheckPointFile(Strand.COMPLEX);
			calc.printSummary( getOputputFilePath() );
		}
		
		else {
			calc.printSummary( getCheckPointFilePath() );
		}
		
		return calc;
	}


	protected PFAbstract deSerializePF( String spName, String path, boolean contSCFlex ) {

		if( !new File(path).exists() ) return null;

		PFAbstract ans = (PFAbstract) ObjectIO.readObject(path, true);
		if(ans != null) {
			name2PF.put(spName, ans);
			// set the panseqsp after de-serializing
			SearchProblem panSeqSP = name2SP.get(getSearchProblemName(contSCFlex, ans.getStrand()));
			ans.setPanSeqSP(panSeqSP);
		}
		return ans;
	}


	protected ArrayList<ArrayList<String>> getSeqsFromFile(String path) {

		ArrayList<ArrayList<String>> ans = new ArrayList<>();
		ArrayList<String> lines = file2List(path);
		if(lines.size() == 0) return ans;

		for( ArrayList<String> seq : strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqList() ) {
			String strSeq = KSAbstract.list1D2String(seq, " ");

			for(String line : lines) {
				if(line.contains(strSeq)) {
					ans.add(seq);
					break;
				}
			}
		}

		return ans;
	}

}
