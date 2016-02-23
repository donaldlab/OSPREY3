package edu.duke.cs.osprey.kstar;

import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;

import edu.duke.cs.osprey.confspace.AllowedSeqs;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.kstar.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.PFAbstract.RunState;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class KSCalc {

	public static enum KSCalcType { SUBSEQUENCE_TREE, PARTIAL_SUBSEQUENCE, FULL_SUBSEQUENCE, FULL_KSTAR }
	protected KSCalcType type = KSCalcType.FULL_KSTAR;

	protected int seqID = -1;
	
	ExpFunction ef = new ExpFunction();

	// each structure has an entry per strand
	HashMap<Integer, PFAbstract> strand2PartitionFunction = new HashMap<>();

	PFAbstract ligand = null;
	PFAbstract protein = null;
	PFAbstract complex = null;

	BigDecimal score = BigDecimal.ZERO;
	int precision = 4;
	static BigDecimal maxValue = new BigDecimal("2e65536");

	public KSCalc( int seqID, 
			HashMap<Integer, AllowedSeqs> allowedSequences,
			HashMap<Integer, SearchProblem> searchProblems,
			HashMap<Integer, PruningControl> pruningControl,
			ConfigFileParser cfp) {

		try {

			double EW = cfp.getParams().getDouble("Ew",0);
			double I0 = cfp.getParams().getBool("imindee",false) ? cfp.getParams().getDouble("Ival",5) : 0;

			this.seqID = seqID;

			HashMap<Integer, ArrayList<String>> strand2Seq = new HashMap<>(); 
			HashMap<Integer, SearchProblem> strand2SearchProblem = new HashMap<>();
			HashMap<Integer, PruningControl> strand2PruningControl = new HashMap<>();

			// set sequence for each state
			strand2Seq.put(Strand.COMPLEX, allowedSequences.get(Strand.COMPLEX).getStrandSeqList().get(seqID));
			strand2Seq.put(Strand.PROTEIN, allowedSequences.get(Strand.PROTEIN).getStrandSeqList().get(seqID));
			strand2Seq.put(Strand.LIGAND, allowedSequences.get(Strand.LIGAND).getStrandSeqList().get(seqID));
			
			ArrayList<String[]> moveableStrands = null;
	        ArrayList<String[]> freeBBZones = null;
	        DEEPerSettings dset = null;

			for( int strand : searchProblems.keySet() ) {
				// set search problem
				ArrayList<String> seq = strand2Seq.get(strand);
				SearchProblem seqSearchProblem = searchProblems.get(strand);

				ArrayList<Integer> pos = new ArrayList<>(); for( int i = 0; i < seq.size(); ++i ) pos.add(i);
				ArrayList<ArrayList<String>> allowedAAs = toListOfLists(seq);
				ArrayList<String> flexibleRes = seqSearchProblem.getFlexibleResiduePositions(seq, pos);
		        
		        moveableStrands = allowedSequences.get(strand).getMoveableStrandTermini();
		        freeBBZones = allowedSequences.get(strand).getFreeBBZoneTermini();
		        dset = allowedSequences.get(strand).getDEEPerSettings();
				
				// new: replace existing search problem
				seqSearchProblem = new SearchProblem( KSImplLinear.getSearchProblemName(strand, seq), 
						cfp.getParams().getValue("PDBNAME"), 
						flexibleRes, allowedAAs, true, true,
						cfp.getParams().getBool("UseEPIC"),
		                new EPICSettings(cfp.getParams()),
		                cfp.getParams().getBool("UseTupExp"),
		                dset, moveableStrands, freeBBZones,
		                cfp.getParams().getBool("useEllipses"),
		                cfp.getParams().getBool("useERef"),
		                cfp.getParams().getBool("AddResEntropy"));
				
				seqSearchProblem.loadEnergyMatrix();

				strand2SearchProblem.put(strand, seqSearchProblem);

				// new
				PruningControl seqPruningControl = cfp.getPruningControl(seqSearchProblem, EW+I0, false, false);
				seqPruningControl.prune();

				strand2PruningControl.put(strand, seqPruningControl);

				// old : set pruning
				//strand2PruningControl.put(strand, pruningControl.get(strand));
			}

			// set partition function for each state
			for( int strand : strand2SearchProblem.keySet() ) {
				setPartitionFunction( strand, 
						PFFactory.getPartitionFunction( PFAbstract.getImplementation(),
								strand2Seq.get(strand), 
								cfp, strand2SearchProblem.get(strand), 
								strand2PruningControl.get(strand), 
								dset, moveableStrands, freeBBZones,
								EW+I0 ) );
			}

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	public KSCalc(int seqID, ArrayList<String> seq, HashMap<Integer, SearchProblem> searchProblems) {
		this.seqID = seqID;
	}


	public KSCalc() {}


	public int getSeqID() {
		return seqID;
	}


	/**
	 * A calculation is valid if it has unpruned residue conformations at all positions
	 * @return
	 */
	public boolean isValid() {

		for( PFAbstract partitionFunction : strand2PartitionFunction.values() ) {
			if( partitionFunction.getNumInitialUnPrunedConfs().equals(BigInteger.ZERO) ) 
				return false;
		}

		return true;
	}


	protected KSCalc expand(KSCalc wtSeq) {

		run(wtSeq);

		return null;
	}


	protected void resetNumMinimizedConfsDuringInterval() {
		// if( ligand != null ) ligand.resetNumMinimizedConfsDuringInterval();
		if( protein != null ) protein.resetNumMinimizedConfsDuringInterval();
		if( complex != null ) complex.resetNumMinimizedConfsDuringInterval();
	}


	protected void run(KSCalc wtSeq) {

		ArrayList<Integer> strands = new ArrayList<>();
		strands.add(Strand.LIGAND);
		strands.add(Strand.PROTEIN);
		strands.add(Strand.COMPLEX);

		for( int strand : strands ) {

			if( getEpsilon() != EApproxReached.FALSE ) return;

			PFAbstract pf = getPartitionFunction(strand);
			ArrayList<String> seq = getSequence(strand);
			
			if( pf.getRunState() == RunState.NOTSTARTED ) {
				System.out.print("\nStarting partition function for "); print(seq, System.out); System.out.println();
				pf.start();
			}

			if( pf.eAppx == EApproxReached.FALSE ) {

				if(PFAbstract.getInterval() == PFAbstract.getMaxInterval()) { 
					System.out.print("\nComputing partition function for "); print(seq, System.out); System.out.println();
					pf.runToCompletion();
				}

				else {
					System.out.print("\nResuming partition function for "); print(seq, System.out); System.out.println();
					pf.resume();
				}
			}
			
			if( wtSeq != null && strand != Strand.COMPLEX ) {
				if( !unboundIsStableWRT(wtSeq, strand) ) {
					
					pf.eAppx = EApproxReached.NOT_STABLE;
					
					return;
				}
			}
		}
	}


	public BigDecimal getKStarScore() {

		BigDecimal divisor = protein.getQStar().multiply( ligand.getQStar() );

		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) score = maxValue;

		else score = complex.getQStar().divide( divisor, precision );

		return score;
	}


	public PruningControl getPruningControl( int strand ) {
		PFAbstract pf = getPartitionFunction(strand);
		return pf == null ? null : pf.pc;
	}


	public PFAbstract getPartitionFunction( int strand ) {
		PFAbstract pf = strand2PartitionFunction.get( strand );
		return pf;
	}


	public void setPartitionFunction( int strand, PFAbstract newPF ) {

		PFAbstract oldPF = getPartitionFunction( strand );
		if( oldPF != null ) oldPF.abort();

		strand2PartitionFunction.put( strand, newPF );

		switch( strand ) {

		case Strand.LIGAND:
			ligand = newPF; break;

		case Strand.PROTEIN:
			protein = newPF; break;

		case Strand.COMPLEX:
			complex = newPF; break;

		default:
			break;
		}
	}


	protected ArrayList<String> getSequence( int strand ) {
		PFAbstract pf = getPartitionFunction(strand);
		return pf == null ? null : pf.sequence;
	}


	protected String getSequenceSignature( int strand ) {

		ArrayList<String> seq = getSequence(strand);
		if(seq == null) return null;

		return getSequenceSignature(seq);
	}


	public static String getSequenceSignature( ArrayList<String> seq ) {

		StringBuilder ans = new StringBuilder();
		for( int i = 0; i < seq.size(); ++i ) {
			if(i > 0) ans.append(".");
			ans.append(seq.get(i));
		}

		return ans.toString();
	}


	public static ArrayList<ArrayList<String>> toListOfLists(ArrayList<String> list) {

		ArrayList<ArrayList<String>> ans = new ArrayList<>();

		for( String s : list ) {

			ArrayList<String> newList = new ArrayList<>();
			newList.add(s);
			ans.add(newList);
		}

		return ans;
	}


	protected SearchProblem getSearchProblem( int strand ) {
		PFAbstract pf = getPartitionFunction(strand);
		return pf == null ? null : pf.sp;
	}


	public void printSequences(PrintStream out) {
		print(getSequence(Strand.LIGAND), out); out.println();
		print(getSequence(Strand.PROTEIN), out); out.println();
		print(getSequence(Strand.COMPLEX), out); out.println();
	}


	public void printSequence(int strand, PrintStream out) {
		print(getSequence(strand), out);
	}


	public static void print(ArrayList<String> seq, PrintStream out) {
		for(int i = 0; i < seq.size(); ++i) {
			if( i != 0 ) out.print(" ");
			out.print(seq.get(i));
		}
	}


	public BigInteger getNumMinimizedConfs() {
		BigInteger total = BigInteger.ZERO;
		if(complex != null) total = total.add(complex.getNumMinimizedConfs());
		if(protein != null) total = total.add(protein.getNumMinimizedConfs());
		//if(ligand != null) total = total.add(ligand.getNumMinimizedConfs());
		return total;
	}


	public void printSummary(PrintStream out) {

		out.print(getSeqID());

		out.print("\t");
		print(getSequence(Strand.COMPLEX), out);

		out.print("\t");
		out.print( ObjectIO.formatBigDecimal(getKStarScore(), precision) );

		ArrayList<Integer> strands = new ArrayList<>();
		strands.add(Strand.COMPLEX);
		strands.add(Strand.PROTEIN);
		strands.add(Strand.LIGAND);

		BigInteger numConfs = BigInteger.ZERO;
		for( int strand : strands ) {
			PFAbstract pf = getPartitionFunction(strand);
			numConfs = numConfs.add( pf.getNumMinimizedConfs() );
		}
		out.print("\t");
		out.print(numConfs);

		for( int strand : strands ) {
			PFAbstract pf = getPartitionFunction(strand);
			out.print("\t");
			out.print( ObjectIO.formatBigDecimal(pf.getQStar(), precision) );

			out.print("\t");
			out.print(pf.getEffectiveEpsilon());

			out.print("\t");
			out.print(pf.getNumMinimizedConfs());
		}

		out.println();
	}


	public void abort() {
		if( ligand != null && ligand.getRunState() != RunState.NOTSTARTED ) ligand.abort();
		if( protein != null && protein.getRunState() != RunState.NOTSTARTED ) protein.abort();
		if( complex != null && complex.getRunState() != RunState.NOTSTARTED ) complex.abort();
	}


	protected EApproxReached getEpsilon() {

		if( complex == null ) return EApproxReached.FALSE;

		if(ligand.eAppx == EApproxReached.NOT_POSSIBLE 
				|| protein.eAppx == EApproxReached.NOT_POSSIBLE 
				|| complex.eAppx == EApproxReached.NOT_POSSIBLE) 
			return EApproxReached.NOT_POSSIBLE;
		
		if(ligand.eAppx == EApproxReached.NOT_STABLE
				|| protein.eAppx == EApproxReached.NOT_STABLE
				|| complex.eAppx == EApproxReached.NOT_STABLE) 
			return EApproxReached.NOT_STABLE;

		else if(ligand.eAppx == EApproxReached.TRUE
				&& protein.eAppx == EApproxReached.TRUE
				&& complex.eAppx == EApproxReached.TRUE) 
			return EApproxReached.TRUE;

		return EApproxReached.FALSE;
	}


	public BigDecimal getKStarUpperBound() {
		BigDecimal bound;
		BigDecimal divisor = protein.getLowerBound().
				multiply(ligand.getLowerBound());

		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) bound = maxValue;

		else bound = complex.getUpperBound().divide( divisor, precision );

		return bound;
	}


	public double getK_StarLowerBound() {

		// if(complex == null) return maxValue;
		if(complex == null) return ef.log10(maxValue);
		
		BigDecimal pl = complex.getUpperBound();
		
		// if( pl.compareTo(BigDecimal.ZERO) == 0 ) return maxValue;
		if( pl.compareTo(BigDecimal.ZERO) == 0 ) return ef.log10(maxValue);

		BigDecimal p = protein.getLowerBound();

		BigDecimal l = ligand.getLowerBound();

		// BigDecimal ans = ( p.multiply( l ) ).divide( pl, precision );
		// return ans;
		
		return ef.log10(l) + ef.log10(p) - ef.log10(pl);
	}


	public double getEAdjustedK_StarLowerBound() {

		// if(complex == null) return new BigDecimal("2e65536");
		if(complex == null) return ef.log10(maxValue);
		
		BigDecimal pl = complex.getUpperBoundAtEpsilon();
		// if( pl.compareTo(BigDecimal.ZERO) == 0) return maxValue;
		if( pl.compareTo(BigDecimal.ZERO) == 0 ) return ef.log10(maxValue);
		
		BigDecimal p = protein.getLowerBound();

		BigDecimal l = ligand.getLowerBound();

		// BigDecimal ans = ( p.multiply( l ) ).divide( pl, precision );
		// return ans;
		
		return ef.log10(l) + ef.log10(p) - ef.log10(pl);
	}


	public double getScaledK_StarLowerBound() {

		double targetEpsilon = PFAbstract.targetEpsilon;

		double ple = getEffectiveEpsilon(Strand.COMPLEX);
		double pe = getEffectiveEpsilon(Strand.PROTEIN);
		double le = getEffectiveEpsilon(Strand.LIGAND);

		BigDecimal plScale = null; double pleDiff = ple <= targetEpsilon ? 0.0 : 1.0-ple;
		BigDecimal pScale = null; double peDiff = pe <= targetEpsilon ? 0.0 : 1.0-pe;
		BigDecimal lScale = null; double leDiff = le <= targetEpsilon ? 0.0 : 1.0-le;

		if( pleDiff > 0.0 ) plScale = new BigDecimal( (1.0-targetEpsilon)/pleDiff );
		if( peDiff > 0.0 ) pScale = new BigDecimal( (1.0-targetEpsilon)/peDiff );
		if( leDiff > 0.0 ) lScale = new BigDecimal( (1.0-targetEpsilon)/leDiff );

		BigDecimal pl = complex.getUpperBound();
		if( plScale != null ) pl = pl.multiply(plScale);
		// if( pl.compareTo(BigDecimal.ZERO) == 0 ) return maxValue;
		if( pl.compareTo(BigDecimal.ZERO) == 0 ) return ef.log10(maxValue);

		BigDecimal p = protein.getLowerBound();
		if( pScale != null ) p = p.multiply(pScale);

		BigDecimal l = ligand.getLowerBound();
		if( lScale != null ) l = l.multiply(lScale);

		// BigDecimal ans = ( p.multiply( l ) ).divide( pl, precision );
		// return ans;
		
		return ef.log10(l) + ef.log10(p) - ef.log10(pl);
	}


	public double getEffectiveEpsilon( int strand ) {

		switch( strand ) {

		case Strand.COMPLEX:
			return complex.getEffectiveEpsilon();

		case Strand.PROTEIN:
			return protein.getEffectiveEpsilon();

		case Strand.LIGAND:
			return ligand.getEffectiveEpsilon();

		default: 
			throw new RuntimeException("ERROR: invalid strand specified");

		}
	}


	protected BigInteger getNumMinimizedConfsDuringInterval() {

		BigInteger ans = BigInteger.ZERO;

		if( complex != null) ans = ans.add(complex.getNumMinimizedConfsDuringInterval());
		if( protein != null ) ans = ans.add(protein.getNumMinimizedConfsDuringInterval());
		if( ligand != null ) ans = ans.add(ligand.getNumMinimizedConfsDuringInterval());

		resetNumMinimizedConfsDuringInterval();

		return ans;
	}


	protected BigInteger getNumInitialUnPrunedConfs() {

		BigInteger ans = complex.getNumInitialUnPrunedConfs();
		ans = ans.add(protein.getNumInitialUnPrunedConfs());
		ans = ans.add(ligand.getNumInitialUnPrunedConfs());

		return ans;
	}


	protected KSCalcType getType() {
		return type;
	}


	protected boolean isRoot() {
		return ligand != null && protein == null && complex == null;
	}


	protected boolean unboundIsStableWRT(KSCalc wtSeq, int strand) {

		if(this.getPartitionFunction(strand).getUpperBound().
				compareTo( wtSeq.getPartitionFunction(strand).getQStar().
						multiply(PFAbstract.getStabilityThreshold()) ) >= 0)
			return true;

		return false;
	}

}
