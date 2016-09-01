package edu.duke.cs.osprey.kstar;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ConcurrentHashMap;

import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.RunState;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class KSCalc {

	private ConcurrentHashMap<Integer, PFAbstract> strand2PF = null;
	private int seqID = 0;
	private static int precision = 4;

	public KSCalc( int seqID, ConcurrentHashMap<Integer, PFAbstract> pfs ) {
		this.seqID = seqID;
		this.strand2PF = pfs;
	}
	

	public PFAbstract getPF(int strand) {
		return strand2PF.get(strand);
	}


	protected boolean unboundIsStable(PFAbstract wtPF, PFAbstract pf) {

		if(pf.getQStarUpperBound().compareTo( wtPF.getQStar().multiply(PFAbstract.getStabilityThresh()) ) >= 0)
			return true;

		return false;
	}

	
	protected boolean doingKAStar() {
		
		for( PFAbstract pf : strand2PF.values() ) {
			if(!pf.getImpl().equalsIgnoreCase(PFAbstract.getCFGImpl())) return true;
		}
		
		return false;
	}
	

	public EApproxReached getEpsilonStatus() {

		PFAbstract pl = strand2PF.get(KSTermini.COMPLEX);
		PFAbstract p = strand2PF.get(KSTermini.PROTEIN);
		PFAbstract l = strand2PF.get(KSTermini.LIGAND);
		
		if(pl.getEpsilonStatus() == EApproxReached.NOT_POSSIBLE 
				|| p.getEpsilonStatus() == EApproxReached.NOT_POSSIBLE 
				|| l.getEpsilonStatus() == EApproxReached.NOT_POSSIBLE) 
			return EApproxReached.NOT_POSSIBLE;

		if(p.getEpsilonStatus() == EApproxReached.NOT_STABLE 
				|| l.getEpsilonStatus() == EApproxReached.NOT_STABLE)
			return EApproxReached.NOT_STABLE;

		else if(pl.getEpsilonStatus() == EApproxReached.TRUE
				&& p.getEpsilonStatus() == EApproxReached.TRUE
				&& l.getEpsilonStatus() == EApproxReached.TRUE) 
			return EApproxReached.TRUE;

		return EApproxReached.FALSE;
	}

	
	public boolean canContinue() {

		if( !doingKAStar() ) {
			if( getEpsilonStatus() == EApproxReached.NOT_POSSIBLE || getEpsilonStatus() == EApproxReached.NOT_STABLE )
				return false;
			
			return true;
		}

		// doing kastar
		if( getPF(KSTermini.COMPLEX).getEpsilonStatus() == EApproxReached.NOT_POSSIBLE )
			return false;
		
		return true;
	}
	
	/*
	public boolean canContinue() {
		
		if( getEpsilonStatus() == EApproxReached.NOT_POSSIBLE 
				|| getEpsilonStatus() == EApproxReached.NOT_STABLE )
			return false;
		
		return true;
	}
	*/
	

	public void runPF(PFAbstract pf, PFAbstract wtPF, boolean complete, boolean stabilityCheck) {
		
		if( pf.getEpsilonStatus() != EApproxReached.FALSE ) return;
		
		// this method shoudl only be called directly for non-complex strands
		if( pf.getRunState() == RunState.NOTSTARTED ) {
			System.out.println("\n" + pf.getImpl() + ": Initializing partition function for " + KSAbstract.list1D2String(pf.getSequence(), " ") + " " + pf.getFlexibility());
			pf.start();
		}

		if( pf.getEpsilonStatus() == EApproxReached.FALSE ) {
			if(complete) {
				System.out.println("\n" + pf.getImpl() + ": Computing partition function for " + KSAbstract.list1D2String(pf.getSequence(), " ") + " " + pf.getFlexibility());
				pf.runToCompletion();
			}

			else {
				System.out.println("\n" + pf.getImpl() + ": Resuming partition function for " + KSAbstract.list1D2String(pf.getSequence(), " ") + " " + pf.getFlexibility());
				pf.runSlice(KSAbstract.checkpointInterval);
			}
		}

		if( pf.getEpsilonStatus() != EApproxReached.FALSE )
			System.out.println("\n" + pf.getImpl() + ": Completed partition function for " + KSAbstract.list1D2String(pf.getSequence(), " ")  + " " + pf.getFlexibility() + "\n");
		
		if( getEpsilonStatus() == EApproxReached.NOT_POSSIBLE ) return;

		if( stabilityCheck && !unboundIsStable(wtPF, pf) ) {
			pf.setEpsilonStatus( EApproxReached.NOT_STABLE );
			System.out.println("\nSequence " + KSAbstract.list1D2String(pf.getSequence(), " ") + " " + pf.getFlexibility() + " is unstable\n");
			return;
		}

	}


	public void run(KSCalc wtKSCalc, boolean forceRun, boolean stabilityCheck) {

		ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(KSTermini.LIGAND, 
				KSTermini.PROTEIN, KSTermini.COMPLEX));

		for( int strand : strands ) {
			if( !forceRun && getEpsilonStatus() != EApproxReached.FALSE ) return;

			boolean complete = KSAbstract.doCheckPoint && strand == KSTermini.COMPLEX ? false : true;
			PFAbstract wtPF = wtKSCalc == null ? null : wtKSCalc.getPF(strand);
			if(wtPF == null || strand == KSTermini.COMPLEX) stabilityCheck = false;
			runPF(getPF(strand), wtPF, complete, stabilityCheck);
		}
	}


	public void serializePFs() {
		for(int strand : strand2PF.keySet()) {
			serializePF(strand);
		}
	}


	public void serializePF( int strand ) {

		PFAbstract pf = strand2PF.get(strand);

		if(pf.getRunState() == RunState.NOTSTARTED) return;

		System.out.print("\nSerializing " + pf.getCheckPointPath() + " ... " );

		pf.abort(false);
		
		if(KSAbstract.doCheckPoint) pf.setPanSeqSP(null);
		
		ObjectIO.writeObject(pf, pf.getCheckPointPath());

		if( pf.getEpsilonStatus() == EApproxReached.FALSE ) {
			// pf state has not been cleaned up
			pf.writeTopConfs();
		}

		System.out.println("Done");	
	}


	public void deleteCheckPointFiles() {

		ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(KSTermini.LIGAND, 
				KSTermini.PROTEIN, KSTermini.COMPLEX));

		for( int strand : strands ) {
			deleteCheckPointFile(strand);
		}
	}


	public void deleteCheckPointFile( int strand ) {
		PFAbstract pf = strand2PF.get(strand);
		ObjectIO.delete(pf.getCheckPointPath());
	}

	
	public BigDecimal getKStarScore( boolean useUB ) {

		BigDecimal pl = useUB ? getPF(KSTermini.COMPLEX).getQStarUpperBound() : getPF(KSTermini.COMPLEX).getQStar();
		BigDecimal p = getPF(KSTermini.PROTEIN).getQStar();
		BigDecimal l = getPF(KSTermini.LIGAND).getQStar();
		
		BigDecimal score;
		
		if( doingKAStar() ) {
			// can easily get clashes for rigid rotamers
			if( (pl.multiply(p).multiply(l)).compareTo(BigDecimal.ZERO) == 0 )
				return score = new BigDecimal(Double.MAX_VALUE);
		}
		
		BigDecimal dividend = pl;
		BigDecimal divisor = p.multiply( l );

		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) {
			
			if(dividend.compareTo(BigDecimal.ZERO) != 0) 
				score = new BigDecimal(Double.POSITIVE_INFINITY);
			
			else
				score = BigDecimal.ZERO;
		}

		else score = dividend.divide( divisor, precision );

		return score;
	}

	
	protected double getKStarScoreLog10( boolean useUB ) {

		BigDecimal pl = useUB ? getPF(KSTermini.COMPLEX).getQStarUpperBound() : getPF(KSTermini.COMPLEX).getQStar();
		BigDecimal p = getPF(KSTermini.PROTEIN).getQStar();
		BigDecimal l = getPF(KSTermini.LIGAND).getQStar();

		if( doingKAStar() ) {
			// can easily get clashes for rigid rotamers in p, l
			if( (pl.multiply(p).multiply(l)).compareTo(BigDecimal.ZERO) == 0 )
			//if( (p.multiply(l)).compareTo(BigDecimal.ZERO) == 0 )
				return Double.POSITIVE_INFINITY;
		}
		
		double score = 0.0;

		ExpFunction e = new ExpFunction();

		if( l.compareTo(BigDecimal.ZERO) == 0 && 
				p.compareTo(BigDecimal.ZERO) == 0 && 
				pl.compareTo(BigDecimal.ZERO) == 0 )
			score = 0.0;

		else if( l.compareTo(BigDecimal.ZERO) == 0 || 
				p.compareTo(BigDecimal.ZERO) == 0 ) {
			
			if(pl.compareTo(BigDecimal.ZERO) != 0)
				score = Double.POSITIVE_INFINITY;
			
			else
				score = 0.0;
		}

		else if( pl.compareTo(BigDecimal.ZERO) == 0 )
			score = Double.NEGATIVE_INFINITY;

		else
			score = e.log10(pl) - e.log10(p) - e.log10(l);

		return score;
	}


	private static void printOutputHeader( PrintStream out ) {

		out.print("Seq ID");
		out.print("\t");
		out.print("Sequence");
		out.print("\t");
		out.print("K* Score (Log10)");
		out.print("\t");
		out.print("UB(K*) Score (Log10)");
		out.print("\t");
		out.print("Total # Confs.");
		out.print("\t");
		out.print("Complex Partition Function");
		out.print("\t");
		out.print("Complex Epsilon");
		out.print("\t");
		out.print("Complex # Confs.");
		out.print("\t");
		out.print("Protein Partition Function");
		out.print("\t");
		out.print("Protein Epsilon");
		out.print("\t");
		out.print("Protein # Confs.");
		out.print("\t");
		out.print("Ligand Partition Function");
		out.print("\t");
		out.print("Ligand Epsilon");
		out.print("\t");
		out.print("Ligand # Confs.");
		out.print("\t");
		out.print("# Seqs Created");
		out.print("\t");
		out.print("# Seqs Completed");
		out.print("\t");
		out.print("Time (sec)");
		out.println();
	}


	public static void printSummaryHeader( String outFile ) {
		try {

			PrintStream out = new PrintStream(new FileOutputStream(outFile, false));
			printOutputHeader(out);
			out.println();

			out.flush();
			out.close();

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	public void printSummary( String outFile, long startTime, long numCreatedSeqs, long numCompletedSeqs ) {

		try {
			
			PrintStream out = new PrintStream(new FileOutputStream(outFile, true));

			out.print(seqID);

			out.print("\t");
			out.print(KSAbstract.list1D2String(getPF(KSTermini.COMPLEX).getSequence(), " "));

			out.print("\t");
			out.print(getKStarScoreLog10(false));
			
			out.print("\t");
			out.print(getKStarScoreLog10(true));

			ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(KSTermini.COMPLEX, 
					KSTermini.PROTEIN, KSTermini.LIGAND));

			BigInteger numMinConfs = BigInteger.ZERO;
			for( int strand : strands ) {
				PFAbstract pf = getPF(strand);
				numMinConfs = numMinConfs.add( pf.getNumMinimized() );
			}
			out.print("\t");
			out.print(numMinConfs);

			for( int strand : strands ) {
				PFAbstract pf = getPF(strand);
				out.print("\t");
				out.print( ObjectIO.formatBigDecimal(pf.getQStar(), precision) );

				out.print("\t");
				out.print(pf.getEffectiveEpsilon());

				out.print("\t");
				out.print(pf.getNumMinimized4Output());
			}
			
			out.print("\t");
			out.print( numCreatedSeqs );
			
			out.print("\t");
			out.print( numCompletedSeqs );
			
			out.print("\t");
			out.print( (System.currentTimeMillis()-startTime)/1000 );
			
			out.println();

			out.flush();
			out.close();

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	public void deleteSeqFromFile( ArrayList<String> seq, String path ) {

		try {

			String strSeq = KSAbstract.list1D2String(seq, " ");

			ArrayList<String> lines = KSAbstract.file2List(path);

			for( String line : lines ) {
				if( line.contains(strSeq) ) {
					lines.remove(line);
					break;
				}
			}

			// write remaing lines to path
			PrintStream out = new PrintStream(new FileOutputStream(path, false));
			for(String line : lines) out.println(line);
			out.flush(); out.close();

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

}
