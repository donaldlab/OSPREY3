package edu.duke.cs.osprey.kstar;

import java.util.HashMap;

import edu.duke.cs.osprey.confspace.AllowedSeqs;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.ConfigFileParser;

public class KSImplLinear2 extends KSAbstract {

	public KSImplLinear2(ConfigFileParser cfp) {
		super(cfp);
	}
	
	public void init( HashMap<Integer, AllowedSeqs> strand2AllowedSeqs ) {
		
		this.strand2AllowedSeqs = strand2AllowedSeqs;
		
		// get search problems per strand
		strand2AllSearchProblem.put(Strand.COMPLEX, cfp.getSearchProblem(Strand.COMPLEX, strand2AllowedSeqs.get(Strand.COMPLEX)));
		strand2AllSearchProblem.put(Strand.PROTEIN, cfp.getSearchProblem(Strand.PROTEIN, strand2AllowedSeqs.get(Strand.PROTEIN)));
		strand2AllSearchProblem.put(Strand.LIGAND, cfp.getSearchProblem(Strand.LIGAND, strand2AllowedSeqs.get(Strand.LIGAND)));
		
		printSequences();
		
		createEnergyMatrices();
	}

	@Override
	public String getKSMethod() {
		return "linear";
	}

	@Override
	public void run() {
		
		if(strand2AllowedSeqs == null)
			throw new RuntimeException("ERROR: call init() method on this object before invoking run()");
		
		// run wt seq
		
	}

}
