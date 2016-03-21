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
	public static boolean refinePInterval = false;
	public static boolean doCheckpoint = false;
	protected static long checkpointInterval = 100000;


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


	protected abstract void prepareAllSingleSeqSPs(ArrayList<Boolean> contSCFlexVals);


	public static void setCheckPointInterval( long interval ) {
		if(interval <= 0) interval = 1;
		checkpointInterval = interval;
	}


	public void createEmats(ArrayList<Boolean> contSCFlexVals) {
		// for now, only the pan seqSPs have energy matrices in the file system
		prepareAllPanSeqSPs(contSCFlexVals);

		// for benchmarking, prepare these ahead of time
		// only single sequences are loaded and pruned this time
		if(refinePInterval)
			prepareAllSingleSeqSPs(contSCFlexVals);
	}


	protected void createCheckPointDir() {
		if( !new File(getCheckPointDir()).exists() )
			ObjectIO.makeDir(getCheckPointDir(), false);
	}


	protected void createEmatDir() {
		if( !new File(getEMATdir()).exists() )
			ObjectIO.makeDir(getEMATdir(), cfp.getParams().getBool("deleteematdir", false));
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


	private String getCheckPointName(boolean contSCFlex, int strand) {

		String flexibility = contSCFlex == true ? "min" : "rig";

		return getCheckPointDir() + 
				File.separator + 
				getRunName() + "." + 
				flexibility + "." + 
				Strand.getStrandString(strand);
	}


	public String getSearchProblemName(boolean contSCFlex, int strand, ArrayList<String> seq) {
		return getSearchProblemName(contSCFlex, strand) + "." + arrayList1D2String(seq, ".");
	}


	public String getCheckPointName(boolean contSCFlex, int strand, ArrayList<String> seq) {
		return getCheckPointName(contSCFlex, strand) + "." + arrayList1D2String(seq, ".") + ".checkpoint";
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


	protected void loadAndPruneMatrices() {

		try {

			System.out.println("\nCreating and pruning pan energy matrices\n");
			long begin = System.currentTimeMillis();

			name2SP.keySet().parallelStream().forEach(key -> {

				// for(String key : name2SP.keySet()) {

				SearchProblem sp = name2SP.get(key);

				// single seq matrices created using the fast construction 
				// method are already pruned according to the pruning window
				if(sp.emat == null) {
					sp.loadEnergyMatrix();
					PruningControl pc = cfp.getPruningControl(sp, EW+I0, false, false); 
					pc.prune();
				}

				// for K* we will guard against under-pruning, since this increases
				// the length of our calculation
				if(sp.isSingleSeq() && KSAbstract.refinePInterval) {
					refinePruningInterval(sp);
				}
				//}

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

		if(sp.numUnPruned().compareTo(BigInteger.ZERO) == 0)
			return;

		// if current interval is good enough, return
		double r = sp.numPruned().doubleValue() / sp.numUnPruned().doubleValue();
		if(r >= pRatioLBT && r <= pRatioUBT) return;

		// most amount of pruning; ratio upper bound
		double l = 0.01;
		PruningControl pc = cfp.getPruningControl(sp, l, false, false); pc.prune();
		double lr = sp.numPruned().doubleValue() / sp.numUnPruned().doubleValue();

		// least amount of pruning; ratio lower bound
		double u = 100;
		pc = cfp.getPruningControl(sp, u, false, false); pc.prune();
		double ur = sp.numPruned().doubleValue() / sp.numUnPruned().doubleValue();

		double m = -1, mr = -1;

		// ratio cannot get smaller or bigger, respectively
		if(ur > pRatioUBT || lr < pRatioLBT) {
			pc = cfp.getPruningControl(sp, EW+I0, false, false); pc.prune();
			sp.pruneMat.setPruningInterval(EW+I0);
			return;
		}

		while( Math.abs(l-u) > 0.1 && (lr > pRatioUBT || ur < pRatioLBT)) {

			m = (l+u)/2;
			pc = cfp.getPruningControl(sp, m, false, false); pc.prune();
			mr = sp.numPruned().doubleValue() / sp.numUnPruned().doubleValue();
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
		pc = cfp.getPruningControl(sp, EW+I0, false, false); pc.prune();
		sp.pruneMat.setPruningInterval(EW+I0);
		return;
	}

	
	protected ConcurrentHashMap<Integer, PFAbstract> createPFsForSeq(ArrayList<ArrayList<String>> seqs, 
			ArrayList<Boolean> contSCFlexVals, ArrayList<String> pfImplVals) {

		ConcurrentHashMap<Integer, PFAbstract> ans = new ConcurrentHashMap<>();

		ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(Strand.COMPLEX, Strand.PROTEIN, Strand.LIGAND));
		ArrayList<Integer> indexes = new ArrayList<>();
		for(int i = 0; i < strands.size(); i++) indexes.add(i);

		indexes.parallelStream().forEach((index) -> {
			// for(int index = 0; index < strands.length; ++index) {

			int strand = strands.get(index);
			boolean contSCFlex = contSCFlexVals.get(index);
			String pfImpl = pfImplVals.get(index);

			AllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);

			ArrayList<String> seq = seqs.get(strand);
			String spName = getSearchProblemName(contSCFlex, strand, seq);
			String cpName = getCheckPointName(contSCFlex, strand, seq);

			if( name2PF.get(spName) == null && deSerializePF(spName, cpName) == null ) {

				ArrayList<Integer> flexResIndexes = strandSeqs.getFlexResIndexesFromSeq(seq);

				// create searchproblem
				SearchProblem seqSP = name2SP.get(spName);
				seqSP = seqSP != null ? seqSP : createSingleSeqSPFast(contSCFlex, strand, seq, flexResIndexes);

				// create partition function
				PFAbstract pf = PFFactory.getPartitionFunction( 
						pfImpl,
						seq,
						cpName,
						cfp, 
						seqSP, 
						EW+I0);

				// put partition function in global map
				name2PF.put(seqSP.name, pf);

				// put in local map
				ans.put(strand, pf);

				// get energy matrix
				if(pf.getSearchProblem().emat == null) {
					pf.getSearchProblem().loadEnergyMatrix();

					// get pruning matrix
					pf.getPruningControl(EW+I0).prune();
				}

				if(seqSP.numUnPruned().compareTo(BigInteger.ZERO) == 0) {
					// no conformations in search space, so this cannot give a valid
					// partition function
					pf.setEpsilonStatus(EApproxReached.NOT_POSSIBLE);
				}

				else {	
					if(!prunedSingleSeqs && KSAbstract.refinePInterval)
						refinePruningInterval(seqSP);

					// initialize conf counts for K*
					pf.setNumUnPruned();
					pf.setNumPruned();
				}
			}
			//}
		});

		// get pfs that were already in global map
		for( int index : indexes ) {

			int strand = strands.get(index);
			boolean contSCFlex = contSCFlexVals.get(index);

			ArrayList<String> seq = seqs.get(strand);
			String spName = getSearchProblemName(contSCFlex, strand, seq);

			if(!ans.keySet().contains(strand)) {

				PFAbstract pf = name2PF.get(spName);

				ans.put(strand, pf);
			}
		}

		if(ans.size() != 3)
			throw new RuntimeException("ERROR: returned map must contain three different partition functions");

		return ans;
	}


	protected SearchProblem createSingleSeqSPSlow( boolean contSCFlex, int strand, ArrayList<String> seq ) {

		ArrayList<ArrayList<String>> allowedAAs = arrayList1D2ListOfLists(AllowedSeqs.getAAsFromSeq(seq));
		ArrayList<String> flexibleRes = AllowedSeqs.getFlexResFromSeq(seq);
		ArrayList<String[]> moveableStrands = strand2AllowedSeqs.get(strand).getMoveableStrandTermini();
		ArrayList<String[]> freeBBZones = strand2AllowedSeqs.get(strand).getFreeBBZoneTermini();
		DEEPerSettings dset = strand2AllowedSeqs.get(strand).getDEEPerSettings();

		// create searchproblem
		SearchProblem seqSP = new SearchProblem( 
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

		return seqSP;
	}


	protected SearchProblem createSingleSeqSPFast( boolean contSCFlex, int strand, 
			ArrayList<String> seq, ArrayList<Integer> flexResIndexes ) {

		String panName = getSearchProblemName(contSCFlex, strand);
		SearchProblem panSeqSP = name2SP.get(panName);

		String singleName = getSearchProblemName(contSCFlex, strand, seq);
		SearchProblem seqSP = panSeqSP.singleSeqSearchProblem(singleName, 
				AllowedSeqs.getAAsFromSeq(seq), AllowedSeqs.getFlexResFromSeq(seq), 
				flexResIndexes);

		return seqSP;
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


	protected void prepareAllPanSeqSPs( ArrayList<Boolean> contSCFlexVals ) {

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

		loadAndPruneMatrices();
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
					+ "." + PFAbstract.getImpl()
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


	protected ArrayList<ArrayList<String>> getStrandStringsAtPos(int i) {

		ArrayList<ArrayList<String>> ans = new ArrayList<ArrayList<String>>(Arrays.asList(null, null, null));

		ans.set(Strand.COMPLEX, strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqAtPos(i));
		ans.set(Strand.PROTEIN, strand2AllowedSeqs.get(Strand.PROTEIN).getStrandSeqAtPos(i));
		ans.set(Strand.LIGAND, strand2AllowedSeqs.get(Strand.LIGAND).getStrandSeqAtPos(i));

		ans.trimToSize();
		return ans;
	}


	protected KSCalc computeWTCalc() {
		// compute wt sequence for reference
		ArrayList<ArrayList<String>> strandSeqs = getStrandStringsAtPos(0);		
		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, true, true));
		ArrayList<String> pfImplVals = new ArrayList<String>(Arrays.asList(PFAbstract.getImpl(), PFAbstract.getImpl(), PFAbstract.getImpl()));

		ConcurrentHashMap<Integer, PFAbstract> pfs = createPFsForSeq(strandSeqs, contSCFlexVals, pfImplVals);
		KSCalc calc = new KSCalc(0, pfs);
		calc.run(calc);
		
		//calc.run(calc, KSAbstract.checkpointInterval);
		//calc.serializePFs();
		//calc.printSummary( getCheckPointFilePath(), true );
		//System.exit(1);
			
		if(calc.getEpsilonStatus() != EApproxReached.TRUE)
			throw new RuntimeException("ERROR: could not compute the wild-type sequence to an epsilon value of "
					+ PFAbstract.targetEpsilon + ". Change the value of epsilon." );

		if(doCheckpoint) {

			for( int strand : Arrays.asList(Strand.LIGAND, Strand.PROTEIN, Strand.COMPLEX) ) {
				PFAbstract pf = calc.getPF(strand);
				if(pf.checkPointExists()) return calc;
			}

			calc.serializePFs();

			calc.printSummary( getCheckPointFilePath(), true );
		}

		calc.printSummary( getOputputFilePath(), true );

		return calc;
	}


	protected PFAbstract deSerializePF( String spName, String path ) {

		if( !new File(path).exists() ) return null;

		PFAbstract ans = (PFAbstract) ObjectIO.readObject(path, true);
		if(ans != null) {
			name2PF.put(spName, ans);
		}
		return ans;
	}


	protected ArrayList<ArrayList<String>> getSeqsFromFile(String path) {

		ArrayList<ArrayList<String>> ans = new ArrayList<>();
		ArrayList<String> lines = file2ArrayList(path);
		if(lines.size() == 0) return ans;

		for( ArrayList<String> seq : strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqList() ) {
			String strSeq = KSAbstract.arrayList1D2String(seq, " ");

			for(String line : lines) {
				if(line.contains(strSeq)) {
					ans.add(seq);
					break;
				}
			}
		}


		return ans;
	}


	public static ArrayList<String> file2ArrayList( String path ) {
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

}
