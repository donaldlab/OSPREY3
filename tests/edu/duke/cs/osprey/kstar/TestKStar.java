package edu.duke.cs.osprey.kstar;

import static edu.duke.cs.osprey.TestBase.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import edu.duke.cs.osprey.control.KStarCalculator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.impl.KSImplLinear;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
(??)import junit.framework.TestResult;

@SuppressWarnings("unused")
public class TestKStar {

	KStarCalculator ksc;
	KSConfigFileParser cfp;

	@Test
	public void test2RL0Linear() {
		
		// read config from files
		KSConfigFileParser cfp = new KSConfigFileParser(new String[] {
			"-c",
			"test/2RL0.kstar/cfgKStar.txt", 
			"Dummy command",
			"test/2RL0.kstar/cfgMutSearch.txt",
			"test/2RL0.kstar/cfgSystem.txt"
		});
		cfp.loadData();
		
		// TEMP
		cfp.getParams().setValue("EmatThreads", "2");
		
		// configure parallelism
		ThreadParallelism.setNumThreads(cfp.getParams().getInt("numThreads", ThreadParallelism.getNumThreads()));
		MultiTermEnergyFunction.setNumThreads(ThreadParallelism.getNumThreads());

		// run K*
		KSImplLinear result = (KSImplLinear)new KStarCalculator(cfp).calcKStarScores();
		
		/* TEMP: print values so we can write this test
		printPfuncs(result, KSTermini.LIGAND);
		printPfuncs(result, KSTermini.PROTEIN);
		printPfuncs(result, KSTermini.COMPLEX);
		printScores(result);
		*/
		
		// check the results
		
		// check ligand partition functions
		List<ArrayList<String>> ligandSequences = result.getUniqueSequences(KSTermini.LIGAND);
		assertThat(ligandSequences.size(), is(15));
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(0), "PHE-156 LYS-172 ILE-192 THR-193", "1.4366299573e+30", 0.75706817);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(1), "PHE-156 LYS-172 ILE-192 SER-193", "1.0514592339e+29", 0.93474129);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(2), "PHE-156 LYS-172 ILE-192 ASN-193", "1.3677261447e+29", 0.78862789);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(3), "PHE-156 LYS-172 ALA-192 THR-193", "4.9965807246e+26", 0.70968666);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(4), "PHE-156 LYS-172 VAL-192 THR-193", "1.7220678832e+28", 0.75569656);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(5), "PHE-156 LYS-172 LEU-192 THR-193", "0.0000000000e+00", 1.00000000);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(6), "PHE-156 LYS-172 PHE-192 THR-193", "8.9945381427e+23", 0.80642125);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(7), "PHE-156 LYS-172 TYR-192 THR-193", "2.4279988599e+25", 0.90656133);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(8), "PHE-156 ASP-172 ILE-192 THR-193", "1.8813145890e+20", 0.63691550);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(9), "PHE-156 GLU-172 ILE-192 THR-193", "2.1237717770e+20", 0.67606823);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(10), "TYR-156 LYS-172 ILE-192 THR-193", "1.7589741810e+30", 0.71250374);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(11), "ALA-156 LYS-172 ILE-192 THR-193", "2.4386287844e+28", 0.34252666);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(12), "VAL-156 LYS-172 ILE-192 THR-193", "6.6126604960e+29", 0.35111757);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(13), "ILE-156 LYS-172 ILE-192 THR-193", "1.8451927552e+30", 0.56251952);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(14), "LEU-156 LYS-172 ILE-192 THR-193", "3.5272424941e+27", 0.44815937);
		
		// check protein partition functions
		List<ArrayList<String>> proteinSequences = result.getUniqueSequences(KSTermini.PROTEIN);
		assertThat(proteinSequences.size(), is(11));
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(0), "PHE-649 ASP-650 GLU-651 THR-654", "1.3839744519e+04", 0.83718344);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(1), "PHE-649 ASP-650 GLU-651 SER-654", "3.3892318567e+05", 0.94815709);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(2), "PHE-649 ASP-650 GLU-651 ASN-654", "1.8829515350e+05", 0.93860704);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(3), "PHE-649 ASP-650 GLU-651 GLN-654", "2.8081086470e+05", 0.94783798);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(4), "PHE-649 ASP-650 ASP-651 THR-654", "3.8050232200e+00", 0.82444117);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(5), "PHE-649 GLU-650 GLU-651 THR-654", "5.7365478276e+04", 0.87284129);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(6), "TYR-649 ASP-650 GLU-651 THR-654", "5.4295126031e+03", 0.87368856);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(7), "ALA-649 ASP-650 GLU-651 THR-654", "3.0001898051e+02", 0.61912380);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(8), "VAL-649 ASP-650 GLU-651 THR-654", "6.0481709740e+01", 0.71981628);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(9), "ILE-649 ASP-650 GLU-651 THR-654", "1.5761115327e+02", 0.93157284);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(10), "LEU-649 ASP-650 GLU-651 THR-654", "2.1833023700e+00", 0.76255380);
		
		// check complex partition functions
		List<ArrayList<String>> complexSequences = result.getUniqueSequences(KSTermini.COMPLEX);
		assertThat(complexSequences.size(), is(25));
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(0), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "2.9714327130e+54", 0.94754675);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(1), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 SER-193", "2.5832792724e+54", 0.94900857);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(2), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 ASN-193", "1.4984350067e+53", 0.94334529);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(3), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ALA-192 THR-193", "3.7011844261e+49", 0.94893920);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(4), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 VAL-192 THR-193", "6.9930872998e+51", 0.93915267);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(5), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 LEU-192 THR-193", "0.0000000000e+00", 1.00000000);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(6), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 PHE-192 THR-193", "0.0000000000e+00", 1.00000000);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(7), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 TYR-192 THR-193", "0.0000000000e+00", 1.00000000);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(8), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 ASP-172 ILE-192 THR-193", "2.0442236895e+41", 0.82510709);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(9), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 GLU-172 ILE-192 THR-193", "3.2269951790e+40", 0.91508768);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(10), "PHE-649 ASP-650 GLU-651 THR-654 TYR-156 LYS-172 ILE-192 THR-193", "2.4839784150e+54", 0.94910417);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(11), "PHE-649 ASP-650 GLU-651 THR-654 ALA-156 LYS-172 ILE-192 THR-193", "1.1382689968e+50", 0.94405927);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(12), "PHE-649 ASP-650 GLU-651 THR-654 VAL-156 LYS-172 ILE-192 THR-193", "1.5913688219e+52", 0.94086988);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(13), "PHE-649 ASP-650 GLU-651 THR-654 ILE-156 LYS-172 ILE-192 THR-193", "2.4042353496e+53", 0.94495208);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(14), "PHE-649 ASP-650 GLU-651 THR-654 LEU-156 LYS-172 ILE-192 THR-193", "6.8484233466e+49", 0.94990203);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(15), "PHE-649 ASP-650 GLU-651 SER-654 PHE-156 LYS-172 ILE-192 THR-193", "3.5117578818e+56", 0.94911386);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(16), "PHE-649 ASP-650 GLU-651 ASN-654 PHE-156 LYS-172 ILE-192 THR-193", "2.3786275180e+56", 0.94894420);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(17), "PHE-649 ASP-650 GLU-651 GLN-654 PHE-156 LYS-172 ILE-192 THR-193", "9.0049639822e+56", 0.94846793);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(18), "PHE-649 ASP-650 ASP-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "8.0539282893e+48", 0.92776587);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(19), "PHE-649 GLU-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "6.8261710438e+53", 0.91486228);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(20), "TYR-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "1.0974114934e+50", 0.94680716);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(21), "ALA-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "1.6247493480e+48", 0.91867532);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(22), "VAL-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "5.3988969294e+48", 0.94121579);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(23), "ILE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "6.2552858594e+49", 0.94007129);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(24), "LEU-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "1.0464124107e+47", 0.93700991);
		
		// check K* scores
		assertThat(result.getKStarScoreLog10(0, false), isRelatively(20.174492908832832));
		assertThat(result.getKStarScoreLog10(1, false), isRelatively(21.249250846820104));
		assertThat(result.getKStarScoreLog10(2, false), isRelatively(19.898510688659158));
		assertThat(result.getKStarScoreLog10(3, false), isRelatively(18.728539744695050));
		assertThat(result.getKStarScoreLog10(4, false), isRelatively(19.467490609604866));
		assertThat(result.getKStarScoreLog10(5, false), isRelatively(0.000000000000000));
		assertThat(result.getKStarScoreLog10(6, false), is(Double.NEGATIVE_INFINITY));
		assertThat(result.getKStarScoreLog10(7, false), is(Double.NEGATIVE_INFINITY));
		assertThat(result.getKStarScoreLog10(8, false), isRelatively(16.894938920347900));
		assertThat(result.getKStarScoreLog10(9, false), isRelatively(16.040562398346555));
		assertThat(result.getKStarScoreLog10(10, false), isRelatively(20.008760279772886));
		assertThat(result.getKStarScoreLog10(11, false), isRelatively(17.527971138235486));
		assertThat(result.getKStarScoreLog10(12, false), isRelatively(18.240266546049790));
		assertThat(result.getKStarScoreLog10(13, false), isRelatively(18.973807164410715));
		assertThat(result.getKStarScoreLog10(14, false), isRelatively(18.147027208034224));
		assertThat(result.getKStarScoreLog10(15, false), isRelatively(20.858078367656290));
		assertThat(result.getKStarScoreLog10(16, false), isRelatively(20.944142378630880));
		assertThat(result.getKStarScoreLog10(17, false), isRelatively(21.348723154946050));
		assertThat(result.getKStarScoreLog10(18, false), isRelatively(18.168305528726854));
		assertThat(result.getKStarScoreLog10(19, false), isRelatively(18.918181629120554));
		assertThat(result.getKStarScoreLog10(20, false), isRelatively(16.148263740488176));
		assertThat(result.getKStarScoreLog10(21, false), isRelatively(15.576292722108057));
		assertThat(result.getKStarScoreLog10(22, false), isRelatively(16.793336058493914));
		assertThat(result.getKStarScoreLog10(23, false), isRelatively(17.441315296074393));
		assertThat(result.getKStarScoreLog10(24, false), isRelatively(16.523244077345900));
	}
	
	private void printPfuncs(KSImplLinear result, int strand) {
		List<ArrayList<String>> sequences = result.getUniqueSequences(strand);
		for (int i=0; i<sequences.size(); i++) {
			printPfunc(result.getPartitionFunction(strand, sequences.get(i)), i, strand);
		}
	}
	
	private void printPfunc(PFAbstract pfunc, int sequenceIndex, int strand) {
		String strandName = KSTermini.getTerminiString(strand);
		System.out.println(String.format("checkPfunc(result, KSTermini.%s, %sSequences.get(%d), \"%s\", \"%.10e\", %.8f);",
			strandName.toUpperCase(),
			strandName.toLowerCase(),
			sequenceIndex,
			KSAbstract.list1D2String(pfunc.getSequence(), " "),
			pfunc.getQStar().doubleValue(),
			pfunc.getEffectiveEpsilon()
		));
	}
	
	private void printScores(KSImplLinear result) {
		for (int i=0; i<result.getSequences(KSTermini.COMPLEX).size(); i++) {
			System.out.println(String.format("assertThat(result.getKStarScoreLog10(%d, false), isRelatively(%.15f));", i, result.getKStarScoreLog10(i, false)));
		}
	}
	
	private void checkPfunc(KSImplLinear result, int strand, ArrayList<String> obsSequence, String expSequenceString, String expQStar, double expEpsilon) {
		
		PFAbstract pfunc = result.getPartitionFunction(strand, obsSequence);
		
		assertThat(pfunc.getStrand(), is(strand));
		
		List<String> expSequence = Arrays.asList(expSequenceString.split(" "));
		assertThat(obsSequence, is(expSequence));
		assertThat(pfunc.getSequence(), is(expSequence));
		
		// NOTE: need a fairly large epsilon since the values are based on minimization
		// and the minimization is very sensitive to the initial molecule state
		// which is different depending on whether or not we computed an energy matrix this run, or read it from a file
		assertThat(pfunc.getQStar(), isRelatively(new BigDecimal(expQStar), 1e-6));
		assertThat(pfunc.getEffectiveEpsilon(), isRelatively(expEpsilon, 1e-5));
	}
}
