package edu.duke.cs.osprey.kstar;

import java.net.InetAddress;
import java.util.ArrayList;
import java.util.HashMap;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.ConfigFileParser;

public abstract class KSAbstract {

	ConfigFileParser cfp = null;
	ArrayList<String> wtSeq = null;
	ArrayList<Integer> computedSeqIDs = new ArrayList<>();
	ArrayList<Integer> prunedSeqIDs = new ArrayList<>();
	HashMap<Integer, HashMap<ArrayList<String>, PFAbstract>> strand2PrecomputedPFs = new HashMap<>();
	
	
	public KSAbstract( ConfigFileParser cfp ) {
		this.cfp = cfp;
		
		strand2PrecomputedPFs.put(Strand.COMPLEX, new HashMap<ArrayList<String>, PFAbstract>());
		strand2PrecomputedPFs.put(Strand.PROTEIN, new HashMap<ArrayList<String>, PFAbstract>());
		strand2PrecomputedPFs.put(Strand.LIGAND, new HashMap<ArrayList<String>, PFAbstract>());
	}
	
	protected abstract String getMethod();
	
	protected String getOputputFileName(String in) {
		
		String ans = null;
		
		try {
			
			ans = in 
					+ "." + InetAddress.getLocalHost().getHostName()
					+ "." + getMethod()
					+ "." + cfp.getParams().getValue("pFuncMethod", "1npcpmcache")
					+ ".txt";
			
		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
		
		return ans;
	}
	
	protected boolean getConcurrentNodeExpansion() {
		
		String expansionType = cfp.getParams().getValue("KStarExpansion", "serial");
		
		return expansionType.equalsIgnoreCase("serial") ? false : true;
	}
	
	protected String getEMATdir() {
		
		return cfp.getParams().getValue("ematdir", "emat");
		
	}
	
	protected ArrayList<String> getWTSeq() {
		if( wtSeq == null ) wtSeq = cfp.getWTSequence();
		return wtSeq;
	}
	
	protected boolean isWT( KSCalc mutSeq ) {
		return mutSeq.complex.getSequence().equals(getWTSeq());
	}
	
	protected void updateMinimizedConfs() {
		
	}
	
	protected void updateTotNumConfs() {
		
	}
}
