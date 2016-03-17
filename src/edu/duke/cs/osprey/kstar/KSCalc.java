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
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class KSCalc {

	private ConcurrentHashMap<Integer, PFAbstract> strand2PF = null;
	private int seqID = 0;
	public static BigDecimal maxValue = new BigDecimal("2e65536");
	private static int precision = 4;

	public KSCalc( int seqID, ConcurrentHashMap<Integer, PFAbstract> pfs ) {
		this.seqID = seqID;
		this.strand2PF = pfs;
	}


	public PFAbstract getPF(int strand) {
		return strand2PF.get(strand);
	}


	protected boolean unboundIsStable(KSCalc wtKSCalc, int strand) {

		PFAbstract wtPf = wtKSCalc.getPF(strand);
		PFAbstract pf = strand2PF.get(strand);
		if(pf.getUpperBound().compareTo( wtPf.getQStar().multiply(PFAbstract.getStabilityThresh()) ) >= 0)
			return true;

		return false;
	}


	public EApproxReached getEpsilonStatus() {

		PFAbstract pl = strand2PF.get(Strand.COMPLEX);
		PFAbstract p = strand2PF.get(Strand.PROTEIN);
		PFAbstract l = strand2PF.get(Strand.LIGAND);

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


	public void run(KSCalc wtKSCalc) {

		ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(Strand.LIGAND, 
				Strand.PROTEIN, Strand.COMPLEX));

		for( int strand : strands ) {

			if( getEpsilonStatus() != EApproxReached.FALSE ) return;

			PFAbstract pf = getPF(strand);
			ArrayList<String> seq = pf.getSequence();

			if( pf.getRunState() == RunState.NOTSTARTED ) {
				System.out.println("\nInitializing partition function for " + KSAbstract.arrayList1D2String(seq, " "));
				pf.start();
			}

			if( pf.getEpsilonStatus() == EApproxReached.FALSE ) {
				System.out.println("\nComputing partition function for " + KSAbstract.arrayList1D2String(seq, " ") + "\n");
				pf.runToCompletion();
			}

			if( getEpsilonStatus() == EApproxReached.NOT_POSSIBLE ) return;

			if( strand != Strand.COMPLEX && !unboundIsStable(wtKSCalc, strand) ) {
				pf.setEpsilonStatus( EApproxReached.NOT_STABLE );
				System.out.println("\nSequence " + KSAbstract.arrayList1D2String(seq, " ") + " is unstable\n");
				return;
			}
		}
	}


	public void run(KSCalc wtKSCalc, long target) {

		ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(Strand.LIGAND, 
				Strand.PROTEIN, Strand.COMPLEX));

		for( int strand : strands ) {

			if( getEpsilonStatus() != EApproxReached.FALSE ) return;

			PFAbstract pf = getPF(strand);
			ArrayList<String> seq = pf.getSequence();

			if( pf.getRunState() == RunState.NOTSTARTED ) {
				System.out.println("\nInitializing partition function for " + KSAbstract.arrayList1D2String(seq, " "));
				pf.start();
			}

			if( pf.getEpsilonStatus() == EApproxReached.FALSE ) {

				if(strand != Strand.COMPLEX) { 
					System.out.println("\nComputing partition function for " + KSAbstract.arrayList1D2String(seq, " ") + "\n");
					pf.runToCompletion();
				}

				else {
					System.out.println("\nResuming partition function for " + KSAbstract.arrayList1D2String(seq, " ") + "\n");
					pf.runSlice(target);
				}
			}

			if( getEpsilonStatus() == EApproxReached.NOT_POSSIBLE ) return;

			if( strand != Strand.COMPLEX && !unboundIsStable(wtKSCalc, strand) ) {
				pf.setEpsilonStatus( EApproxReached.NOT_STABLE );
				System.out.println("\nSequence " + KSAbstract.arrayList1D2String(seq, " ") + " is unstable\n");
				return;
			}
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
		ObjectIO.writeObject(pf, pf.getCheckPointPath());

		if( pf.getEpsilonStatus() == EApproxReached.FALSE ) {
			// pf state has not been cleaned up
			pf.writeTopConfs();
		}

		System.out.println("Done");	
	}


	public void deleteCheckPointFiles() {

		ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(Strand.LIGAND, 
				Strand.PROTEIN, Strand.COMPLEX));

		for( int strand : strands ) {
			deleteCheckPointFile(strand);
		}
	}


	public void deleteCheckPointFile( int strand ) {
		PFAbstract pf = strand2PF.get(strand);
		ObjectIO.delete(pf.getCheckPointPath());
	}


	protected BigDecimal getKStarScore() {

		BigDecimal score = BigDecimal.ZERO;

		PFAbstract pl = getPF(Strand.COMPLEX);
		PFAbstract p = getPF(Strand.PROTEIN);
		PFAbstract l = getPF(Strand.LIGAND);

		BigDecimal divisor = p.getQStar().multiply( l.getQStar() );

		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) 
			score = maxValue;

		else score = pl.getQStar().divide( divisor, precision );

		return score;
	}


	private void printOutputHeader( PrintStream out ) {

		out.print("Seq ID");
		out.print("\t");
		out.print("Sequence");
		out.print("\t");
		out.print("K* Score");
		out.print("\t");
		out.print("Total # Confs.");
		out.print("\t");
		out.print("Complex Partition Function");
		out.print("\t");
		out.print("Complex Epsilon");
		out.print("\t");
		out.print("Complex # Minimized Confs.");
		out.print("\t");
		out.print("Protein Partition Function");
		out.print("\t");
		out.print("Protein Epsilon");
		out.print("\t");
		out.print("Protein # Minimized Confs.");
		out.print("\t");
		out.print("Ligand Partition Function");
		out.print("\t");
		out.print("Ligand Epsilon");
		out.print("\t");
		out.print("Ligand # Minimized Confs.");
		out.println();
	}


	public void printSummary( String outFile, boolean printHeader ) {

		try {

			boolean append = printHeader ? false : true;

			PrintStream out = new PrintStream(new FileOutputStream(outFile, append));

			if(printHeader)
				printOutputHeader(out);

			out.print(seqID);

			out.print("\t");
			out.print(KSAbstract.arrayList1D2String(getPF(Strand.COMPLEX).getSequence(), " "));

			out.print("\t");
			out.print( ObjectIO.formatBigDecimal(getKStarScore(), precision) );

			ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(Strand.COMPLEX, 
					Strand.PROTEIN, Strand.LIGAND));

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
				out.print(pf.getNumMinimized());
			}

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

			String strSeq = KSAbstract.arrayList1D2String(seq, " ");

			ArrayList<String> lines = KSAbstract.file2ArrayList(path);

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
