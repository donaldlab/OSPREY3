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

@SuppressWarnings("unused")
public class TestKStar {

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
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(0), "PHE-156 LYS-172 ILE-192 THR-193", "1.4366299409e+30", 0.75706528);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(1), "PHE-156 LYS-172 ILE-192 SER-193", "1.0514591737e+29", 0.93474057);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(2), "PHE-156 LYS-172 ILE-192 ASN-193", "1.3677262102e+29", 0.78862541);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(3), "PHE-156 LYS-172 ALA-192 THR-193", "4.9965807004e+26", 0.70968301);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(4), "PHE-156 LYS-172 VAL-192 THR-193", "1.7220678752e+28", 0.75569340);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(5), "PHE-156 LYS-172 LEU-192 THR-193", "0.0000000000e+00", 1.00000000);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(6), "PHE-156 LYS-172 PHE-192 THR-193", "8.9945381365e+23", 0.80641876);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(7), "PHE-156 LYS-172 TYR-192 THR-193", "2.4279987723e+25", 0.90656018);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(8), "PHE-156 ASP-172 ILE-192 THR-193", "1.8813145868e+20", 0.63691550);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(9), "PHE-156 GLU-172 ILE-192 THR-193", "2.1237717041e+20", 0.67606823);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(10), "TYR-156 LYS-172 ILE-192 THR-193", "1.7589742663e+30", 0.71250024);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(11), "ALA-156 LYS-172 ILE-192 THR-193", "2.4386285171e+28", 0.34251749);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(12), "VAL-156 LYS-172 ILE-192 THR-193", "6.6126604925e+29", 0.35110847);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(13), "ILE-156 LYS-172 ILE-192 THR-193", "1.8451925941e+30", 0.56251382);
		checkPfunc(result, KSTermini.LIGAND, ligandSequences.get(14), "LEU-156 LYS-172 ILE-192 THR-193", "3.5272426952e+27", 0.44815176);
		
		// check protein partition functions
		List<ArrayList<String>> proteinSequences = result.getUniqueSequences(KSTermini.PROTEIN);
		assertThat(proteinSequences.size(), is(11));
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(0), "PHE-649 ASP-650 GLU-651 THR-654", "1.3839746270e+04", 0.83718370);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(1), "PHE-649 ASP-650 GLU-651 SER-654", "3.3892319341e+05", 0.94815717);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(2), "PHE-649 ASP-650 GLU-651 ASN-654", "1.8829515352e+05", 0.93860714);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(3), "PHE-649 ASP-650 GLU-651 GLN-654", "2.8081084728e+05", 0.94783806);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(4), "PHE-649 ASP-650 ASP-651 THR-654", "3.8050232200e+00", 0.82444118);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(5), "PHE-649 GLU-650 GLU-651 THR-654", "5.7365478240e+04", 0.87284151);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(6), "TYR-649 ASP-650 GLU-651 THR-654", "5.4295126690e+03", 0.87368784);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(7), "ALA-649 ASP-650 GLU-651 THR-654", "3.0001898088e+02", 0.61912380);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(8), "VAL-649 ASP-650 GLU-651 THR-654", "6.0481644580e+01", 0.71981649);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(9), "ILE-649 ASP-650 GLU-651 THR-654", "1.5761115359e+02", 0.93157284);
		checkPfunc(result, KSTermini.PROTEIN, proteinSequences.get(10), "LEU-649 ASP-650 GLU-651 THR-654", "2.1833023600e+00", 0.76255380);
		
		// check complex partition functions
		List<ArrayList<String>> complexSequences = result.getUniqueSequences(KSTermini.COMPLEX);
		assertThat(complexSequences.size(), is(25));
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(0), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "2.9714808704e+54", 0.94754644);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(1), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 SER-193", "2.5815043174e+54", 0.94904214);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(2), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 ASN-193", "1.4984223821e+53", 0.94334673);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(3), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ALA-192 THR-193", "3.6992710121e+49", 0.94896455);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(4), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 VAL-192 THR-193", "6.9930812772e+51", 0.93915343);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(5), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 LEU-192 THR-193", "0.0000000000e+00", 1.00000000);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(6), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 PHE-192 THR-193", "0.0000000000e+00", 1.00000000);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(7), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 TYR-192 THR-193", "0.0000000000e+00", 1.00000000);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(8), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 ASP-172 ILE-192 THR-193", "2.0442232676e+41", 0.82510719);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(9), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 GLU-172 ILE-192 THR-193", "3.2269951524e+40", 0.91508772);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(10), "PHE-649 ASP-650 GLU-651 THR-654 TYR-156 LYS-172 ILE-192 THR-193", "2.4402030733e+54", 0.94622697);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(11), "PHE-649 ASP-650 GLU-651 THR-654 ALA-156 LYS-172 ILE-192 THR-193", "1.1382664631e+50", 0.94405991);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(12), "PHE-649 ASP-650 GLU-651 THR-654 VAL-156 LYS-172 ILE-192 THR-193", "1.5903923683e+52", 0.94090457);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(13), "PHE-649 ASP-650 GLU-651 THR-654 ILE-156 LYS-172 ILE-192 THR-193", "2.4003079715e+53", 0.94503717);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(14), "PHE-649 ASP-650 GLU-651 THR-654 LEU-156 LYS-172 ILE-192 THR-193", "6.8484619549e+49", 0.94990228);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(15), "PHE-649 ASP-650 GLU-651 SER-654 PHE-156 LYS-172 ILE-192 THR-193", "3.5049487131e+56", 0.94920794);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(16), "PHE-649 ASP-650 GLU-651 ASN-654 PHE-156 LYS-172 ILE-192 THR-193", "2.3835564150e+56", 0.94884425);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(17), "PHE-649 ASP-650 GLU-651 GLN-654 PHE-156 LYS-172 ILE-192 THR-193", "9.0336655117e+56", 0.94831250);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(18), "PHE-649 ASP-650 ASP-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "8.0539291616e+48", 0.92776785);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(19), "PHE-649 GLU-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "6.8261719992e+53", 0.91486638);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(20), "TYR-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "1.0974127664e+50", 0.94680725);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(21), "ALA-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "1.6247493349e+48", 0.91867560);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(22), "VAL-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "5.5252601391e+48", 0.93992280);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(23), "ILE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "6.1974332024e+49", 0.94059299);
		checkPfunc(result, KSTermini.COMPLEX, complexSequences.get(24), "LEU-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "1.0201459863e+47", 0.93849399);
		
		// check K* scores
		assertThat(result.getKStarScoreLog10(0, false), isRelatively(20.174499897285763));
		assertThat(result.getKStarScoreLog10(1, false), isRelatively(21.248952313159236));
		assertThat(result.getKStarScoreLog10(2, false), isRelatively(19.898506953918336));
		assertThat(result.getKStarScoreLog10(3, false), isRelatively(18.728315115076917));
		assertThat(result.getKStarScoreLog10(4, false), isRelatively(19.467490182640140));
		assertThat(result.getKStarScoreLog10(5, false), isRelatively(0.000000000000000));
		assertThat(result.getKStarScoreLog10(6, false), is(Double.NEGATIVE_INFINITY));
		assertThat(result.getKStarScoreLog10(7, false), is(Double.NEGATIVE_INFINITY));
		assertThat(result.getKStarScoreLog10(8, false), isRelatively(16.894938776280530));
		assertThat(result.getKStarScoreLog10(9, false), isRelatively(16.040562354710698));
		assertThat(result.getKStarScoreLog10(10, false), isRelatively(20.001038355879018));
		assertThat(result.getKStarScoreLog10(11, false), isRelatively(17.527970164203865));
		assertThat(result.getKStarScoreLog10(12, false), isRelatively(18.239999929263895));
		assertThat(result.getKStarScoreLog10(13, false), isRelatively(18.973097136499895));
		assertThat(result.getKStarScoreLog10(14, false), isRelatively(18.147029576671514));
		assertThat(result.getKStarScoreLog10(15, false), isRelatively(20.857235464308530));
		assertThat(result.getKStarScoreLog10(16, false), isRelatively(20.945041380135720));
		assertThat(result.getKStarScoreLog10(17, false), isRelatively(21.350105212711505));
		assertThat(result.getKStarScoreLog10(18, false), isRelatively(18.168305580709905));
		assertThat(result.getKStarScoreLog10(19, false), isRelatively(18.918181695122435));
		assertThat(result.getKStarScoreLog10(20, false), isRelatively(16.148264243918618));
		assertThat(result.getKStarScoreLog10(21, false), isRelatively(15.576292723005167));
		assertThat(result.getKStarScoreLog10(22, false), isRelatively(16.803384225115607));
		assertThat(result.getKStarScoreLog10(23, false), isRelatively(17.437279993418002));
		assertThat(result.getKStarScoreLog10(24, false), isRelatively(16.512203527644587));
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
		
		assertThat(pfunc.getQStar(), isRelatively(new BigDecimal(expQStar), 1e-10));
		assertThat(pfunc.getEffectiveEpsilon(), isRelatively(expEpsilon));
	}
}
