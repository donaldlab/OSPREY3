package edu.duke.cs.osprey.kstar;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.kstar.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.PFAbstract.RunState;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class KSCalc {

	private HashMap<Integer, PFAbstract> strand2PF = null;
	private int seqID = 0;
	private boolean headerPrinted = false;
	public static BigDecimal maxValue = new BigDecimal("2e65536");
	private static int precision = 4;

	public KSCalc( int seqID, HashMap<Integer, PFAbstract> strand2PF ) {
		this.seqID = seqID;
		this.strand2PF = strand2PF;
	}

	protected PFAbstract getPF(int strand) {
		return strand2PF.get(strand);
	}

	protected boolean unboundIsStable(KSCalc wtKSCalc, int strand) {

		PFAbstract wtPf = wtKSCalc.getPF(strand);
		PFAbstract pf = strand2PF.get(strand);
		if(pf.getUpperBound().compareTo( wtPf.getQStar().multiply(PFAbstract.getStabilityThreshold()) ) >= 0)
			return true;

		return false;
	}

	protected EApproxReached getEpsilonStatus() {

		PFAbstract pl = strand2PF.get(Strand.COMPLEX);
		PFAbstract p = strand2PF.get(Strand.PROTEIN);
		PFAbstract l = strand2PF.get(Strand.LIGAND);

		if(pl.eAppx == EApproxReached.NOT_POSSIBLE 
				|| p.eAppx == EApproxReached.NOT_POSSIBLE 
				|| l.eAppx == EApproxReached.NOT_POSSIBLE) 
			return EApproxReached.NOT_POSSIBLE;

		if(p.eAppx == EApproxReached.NOT_STABLE 
				|| l.eAppx == EApproxReached.NOT_STABLE)
			return EApproxReached.NOT_STABLE;

		else if(pl.eAppx == EApproxReached.TRUE
				&& p.eAppx == EApproxReached.TRUE
				&& l.eAppx == EApproxReached.TRUE) 
			return EApproxReached.TRUE;

		return EApproxReached.FALSE;
	}

	protected void run(KSCalc wtKSCalc) {

		ArrayList<Integer> strands = new ArrayList<>();
		strands.add(Strand.PROTEIN);
		strands.add(Strand.LIGAND);
		strands.add(Strand.COMPLEX);

		for( int strand : strands ) {

			if( getEpsilonStatus() != EApproxReached.FALSE ) 
				return;

			PFAbstract pf = getPF(strand);
			ArrayList<String> seq = pf.getSequence();

			if( pf.getRunState() == RunState.NOTSTARTED ) {

				System.out.println("\nInitializing partition function for " + KSAbstract.arrayList1D2String(seq, " "));
				pf.start();
			}

			if( pf.eAppx == EApproxReached.FALSE ) {

				if(PFAbstract.getInterval() == PFAbstract.getMaxInterval()) { 

					System.out.println("\nComputing partition function for " + KSAbstract.arrayList1D2String(seq, " "));
					pf.runToCompletion();
				}

				else {

					System.out.println("\nResuming partition function for " + KSAbstract.arrayList1D2String(seq, " "));
					pf.resume();
				}
			}

			if( strand != Strand.COMPLEX && !unboundIsStable(wtKSCalc, strand) ) {

				pf.eAppx = EApproxReached.NOT_STABLE;
				return;
			}
		}

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

		headerPrinted = true;
	}

	public void printSummary( String outFile ) {

		try {

			PrintStream out = new PrintStream(new FileOutputStream(outFile, true));

			if(!headerPrinted) 
				printOutputHeader(out);

			out.print(seqID);

			out.print("\t");
			KSAbstract.arrayList1D2String(getPF(Strand.COMPLEX).getSequence(), " ");

			out.print("\t");
			out.print( ObjectIO.formatBigDecimal(getKStarScore(), precision) );

			ArrayList<Integer> strands = new ArrayList<>();
			strands.add(Strand.COMPLEX);
			strands.add(Strand.PROTEIN);
			strands.add(Strand.LIGAND);

			BigInteger numMinConfs = BigInteger.ZERO;
			for( int strand : strands ) {
				PFAbstract pf = getPF(strand);
				numMinConfs = numMinConfs.add( pf.getNumMinimizedConfs() );
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
				out.print(pf.getNumMinimizedConfs());
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

}
