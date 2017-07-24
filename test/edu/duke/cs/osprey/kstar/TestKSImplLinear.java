package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.KStarCalculator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.impl.KSImplLinear;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;

public class TestKSImplLinear {

	private void testLinear(KSConfigFileParser cfp) {
		
		double targetEpsilon = cfp.params.getDouble("epsilon");
		
		// run K*
		KSImplLinear result = (KSImplLinear)new KStarCalculator(cfp).calcKStarScores();
		
		/* TEMP: print values so we can write this test
		printPfuncs(result, KS1);
		printPfuncs(result, KS0);
		printPfuncs(result, KS2);
		*/
		
		// check the results
		// NOTE: expected values were calculated with targetEpsilon=0.05
		
		// check ligand partition functions
		List<ArrayList<String>> ligandSequences = result.getUniqueSequences(1);
		assertThat(ligandSequences.size(), is(15));
		checkPfunc(result, 1, ligandSequences.get(0), "PHE-156 LYS-172 ILE-192 THR-193", "4.3486162283e+30", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(1), "PHE-156 LYS-172 ILE-192 SER-193", "1.1137761298e+30", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(2), "PHE-156 LYS-172 ILE-192 ASN-193", "4.6813161378e+29", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(3), "PHE-156 LYS-172 ALA-192 THR-193", "1.5510745584e+27", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(4), "PHE-156 LYS-172 VAL-192 THR-193", "5.7469787121e+28", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(5), "PHE-156 LYS-172 LEU-192 THR-193", EApproxReached.NOT_POSSIBLE);
		checkPfunc(result, 1, ligandSequences.get(6), "PHE-156 LYS-172 PHE-192 THR-193", "2.8245771923e+24", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(7), "PHE-156 LYS-172 TYR-192 THR-193", "1.4301116312e+26", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(8), "PHE-156 ASP-172 ILE-192 THR-193", "4.2890115335e+20", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(9), "PHE-156 GLU-172 ILE-192 THR-193", "4.8319041432e+20", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(10), "TYR-156 LYS-172 ILE-192 THR-193", "4.5794607486e+30", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(11), "ALA-156 LYS-172 ILE-192 THR-193", "3.2901636943e+28", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(12), "VAL-156 LYS-172 ILE-192 THR-193", "9.0045148982e+29", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(13), "ILE-156 LYS-172 ILE-192 THR-193", "3.3840619669e+30", targetEpsilon);
		checkPfunc(result, 1, ligandSequences.get(14), "LEU-156 LYS-172 ILE-192 THR-193", "5.3110034153e+27", targetEpsilon);

		
		// check protein partition functions
		List<ArrayList<String>> proteinSequences = result.getUniqueSequences(0);
		assertThat(proteinSequences.size(), is(11));
		checkPfunc(result, 0, proteinSequences.get(0), "PHE-649 ASP-650 GLU-651 THR-654", "4.3053687181e+04", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(1), "PHE-649 ASP-650 GLU-651 SER-654", "3.2226637016e+06", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(2), "PHE-649 ASP-650 GLU-651 ASN-654", "1.5125751083e+06", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(3), "PHE-649 ASP-650 GLU-651 GLN-654", "2.5770817910e+06", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(4), "PHE-649 ASP-650 ASP-651 THR-654", "1.2061571830e+01", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(5), "PHE-649 GLU-650 GLU-651 THR-654", "1.9981960274e+05", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(6), "TYR-649 ASP-650 GLU-651 THR-654", "1.6852241686e+04", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(7), "ALA-649 ASP-650 GLU-651 THR-654", "6.1007791897e+02", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(8), "VAL-649 ASP-650 GLU-651 THR-654", "1.2701407840e+02", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(9), "ILE-649 ASP-650 GLU-651 THR-654", "5.9845437720e+02", targetEpsilon);
		checkPfunc(result, 0, proteinSequences.get(10), "LEU-649 ASP-650 GLU-651 THR-654", "4.5709691200e+00", targetEpsilon);

		
		// check complex partition functions
		List<ArrayList<String>> complexSequences = result.getUniqueSequences(2);
		assertThat(complexSequences.size(), is(25));
		checkPfunc(result, 2, complexSequences.get(0), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "3.5517551983e+54", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(1), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 SER-193", "3.1725379840e+54", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(2), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 ASN-193", "1.6135714647e+53", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(3), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ALA-192 THR-193", "5.6751939710e+49", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(4), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 VAL-192 THR-193", "8.2838823573e+51", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(5), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 LEU-192 THR-193", EApproxReached.FALSE);
		checkPfunc(result, 2, complexSequences.get(6), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 PHE-192 THR-193", EApproxReached.NOT_POSSIBLE);
		checkPfunc(result, 2, complexSequences.get(7), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 TYR-192 THR-193", EApproxReached.NOT_POSSIBLE);
		checkPfunc(result, 2, complexSequences.get(8), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 ASP-172 ILE-192 THR-193", "2.2004951030e+41", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(9), "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 GLU-172 ILE-192 THR-193", "3.5518926769e+40", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(10), "PHE-649 ASP-650 GLU-651 THR-654 TYR-156 LYS-172 ILE-192 THR-193", "2.6478756391e+54", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(11), "PHE-649 ASP-650 GLU-651 THR-654 ALA-156 LYS-172 ILE-192 THR-193", "1.3859773956e+50", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(12), "PHE-649 ASP-650 GLU-651 THR-654 VAL-156 LYS-172 ILE-192 THR-193", "1.8764812326e+52", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(13), "PHE-649 ASP-650 GLU-651 THR-654 ILE-156 LYS-172 ILE-192 THR-193", "2.6145750849e+53", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(14), "PHE-649 ASP-650 GLU-651 THR-654 LEU-156 LYS-172 ILE-192 THR-193", "8.0098311188e+49", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(15), "PHE-649 ASP-650 GLU-651 SER-654 PHE-156 LYS-172 ILE-192 THR-193", "4.3105596294e+56", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(16), "PHE-649 ASP-650 GLU-651 ASN-654 PHE-156 LYS-172 ILE-192 THR-193", "2.8654652076e+56", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(17), "PHE-649 ASP-650 GLU-651 GLN-654 PHE-156 LYS-172 ILE-192 THR-193", "1.0786109248e+57", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(18), "PHE-649 ASP-650 ASP-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "1.1091956414e+49", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(19), "PHE-649 GLU-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "8.0725173333e+53", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(20), "TYR-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "2.2258393553e+50", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(21), "ALA-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "2.7792831743e+48", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(22), "VAL-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "9.1350184819e+48", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(23), "ILE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "8.0042346894e+49", targetEpsilon);
		checkPfunc(result, 2, complexSequences.get(24), "LEU-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "1.8064201386e+47", targetEpsilon);
	}
	
	@SuppressWarnings("unused")
	private void printPfuncs(KSImplLinear result, int strand) {
		List<ArrayList<String>> sequences = result.getUniqueSequences(strand);
		for (int i=0; i<sequences.size(); i++) {
			printPfunc(result.getPartitionFunction(strand, sequences.get(i)), i, strand);
		}
	}
	
	private void printPfunc(PFAbstract pfunc, int sequenceIndex, int strand) {
		String strandName = "Strand"+strand;
		System.out.print(String.format("checkPfunc(result, KSTermini.%s, %sSequences.get(%d), \"%s\"",
			strandName.toUpperCase(),
			strandName.toLowerCase(),
			sequenceIndex,
			KSAbstract.list1D2String(pfunc.getSequence(), " ")
		));
		if (pfunc.getEpsilonStatus() == EApproxReached.TRUE) {
			System.out.print(String.format(", \"%.10e\", targetEpsilon", pfunc.getQStar().doubleValue()));
		} else {
			System.out.print(", EApproxReached." + pfunc.getEpsilonStatus().name());
		}
		System.out.println(");");
	}
	
	private void checkPfunc(KSImplLinear result, int strand, ArrayList<String> obsSequence, String expSequenceString, String expQStar, double targetEpsilon) {
		
		checkPfunc(result, strand, obsSequence, expSequenceString, EApproxReached.TRUE);
		
		PFAbstract pfunc = result.getPartitionFunction(strand, obsSequence);
		
		assertThat(pfunc.getEpsilonStatus(), is(EApproxReached.TRUE));
		
		// we should be within epsilon of the higher precision value for q*
		double minQStar = new BigDecimal(expQStar).doubleValue()*(1.0 - targetEpsilon);
		assertThat(pfunc.getQStar().doubleValue(), greaterThanOrEqualTo(minQStar));
		
		assertThat(pfunc.getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
	}
	
	private void checkPfunc(KSImplLinear result, int strand, ArrayList<String> obsSequence, String expSequenceString, EApproxReached epsilonStatus) {
		
		PFAbstract pfunc = result.getPartitionFunction(strand, obsSequence);
		List<String> expSequence = Arrays.asList(expSequenceString.split(" "));
		
		assertThat(pfunc.getStrand(), is(strand));
		assertThat(obsSequence, is(expSequence));
		assertThat(pfunc.getSequence(), is(expSequence));
		assertThat(pfunc.getEpsilonStatus(), is(epsilonStatus));
	}
	
	private KSConfigFileParser make2RL0Config() {
		
		// read config from files
		KSConfigFileParser cfp = new KSConfigFileParser(ConfigFileParser.makeFromFilePaths(
			"examples/2RL0.kstar/cfgMutSearch.txt",
			"examples/2RL0.kstar/cfgSystem.txt"
		));
		cfp.loadData();
		
		// override file-based config
		
		// I'm guessing most people have at least two cores, so compute the energy matrix a bit faster
		cfp.params.setValue("EmatThreads", "2");
		
		// this test takes several minutes at the config file's value of e=0.95,
		// but it only takes about two minutes at e=0.99
		cfp.params.setValue("epsilon", "0.99");
		
		return cfp;
	}
	
	@Test
	public void test2RL0LinearParallel0() {
		
		// configure parallel0 parallelism
		ThreadParallelism.setNumThreadsIfPossible(1);
		MultiTermEnergyFunction.setNumThreads(1);
		
		KSConfigFileParser cfp = make2RL0Config();
		testLinear(cfp);
	}
	
	@Test
	public void test2RL0LinearParallelConf() {
		KSConfigFileParser cfp = make2RL0Config();
		cfp.params.setValue("kStarPFuncMethod", "parallelConf");
		testLinear(cfp);
	}
	
	@Test
	public void test2RL0LinearParallelConf2Threads() {
		KSConfigFileParser cfp = make2RL0Config();
		cfp.params.setValue("kStarPFuncMethod", "parallelConf");
		cfp.params.setValue("MinimizationThreads", "2");
		testLinear(cfp);
	}
	
	@Test
	public void test2RL0LinearParallelConf1Gpu() {
		KSConfigFileParser cfp = make2RL0Config();
		cfp.params.setValue("kStarPFuncMethod", "parallelConf");
		cfp.params.setValue("MinimizationGpus", "1");
		testLinear(cfp);
	}
	
	@Test
	public void test2RL0LinearParallelConf1Gpu4Streams() {
		KSConfigFileParser cfp = make2RL0Config();
		cfp.params.setValue("kStarPFuncMethod", "parallelConf");
		cfp.params.setValue("MinimizationGpus", "1");
		cfp.params.setValue("MinimizationStreamsPerGpu", "4");
		testLinear(cfp);
	}
	
	@Test
	public void test2RL0LinearParallelConf1Gpu16Streams() {
		KSConfigFileParser cfp = make2RL0Config();
		cfp.params.setValue("kStarPFuncMethod", "parallelConf");
		cfp.params.setValue("MinimizationGpus", "1");
		cfp.params.setValue("MinimizationStreamsPerGpu", "16");
		testLinear(cfp);
	}
}
