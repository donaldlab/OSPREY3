package edu.duke.cs.osprey.ewakstar;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayList;

import edu.duke.cs.osprey.ewakstar.PartitionFunction;
import edu.duke.cs.osprey.multistatekstar.KStarScore;
import edu.duke.cs.osprey.tools.BigDecimalUtil;

/**
*
* @author Anna Lowegard(anna.lowegard@duke.edu)
* Adegoke Ojewole (ao68@duke.edu)
*/

public class EWAKScore {

	private ArrayList<PartitionFunction> pfs;
	private BigDecimal score;
	
	public EWAKScore() {
		pfs = new ArrayList<>();
		score = null;
	}
	
	public void add(PartitionFunction pf) {
		pfs.add(pf);
		pfs.trimToSize();
	}
	
	public BigDecimal getDenominator() {
		PartitionFunction pf;
		BigDecimal ans = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		for(int strand=0;strand<pfs.size()-1;++strand) {
			pf = pfs.get(strand);
			if(pf==null || pf.getZ().compareTo(BigDecimal.ZERO)==0) {
				return BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			}
			ans = ans.multiply(pf.getZ());
		}
		return ans;
	}
	
	public BigDecimal getNumerator() {
		BigDecimal ans = pfs.get(pfs.size()-1).getZ();
		if(ans == null) ans = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
		return ans;
	}
	
	public BigDecimal toLog10(BigDecimal val) {
		if(val.compareTo(BigDecimal.ZERO)==0)
			return KStarScore.MIN_VALUE;

		return BigDecimalUtil.log10(val);
	}
	
	public BigDecimal getScoreLog10() {
		if(score != null) {
			return score;
		}
		
		BigDecimal den = getDenominator();
		if(den.compareTo(BigDecimal.ZERO) == 0) return toLog10(den);
		BigDecimal ans = getNumerator().setScale(64, RoundingMode.HALF_UP).divide(den, RoundingMode.HALF_UP);
		score = toLog10(ans);
		return score;
	}
	
	public BigDecimal getScore() {
		if(score != null) {
			return score;
		}
		
		BigDecimal den = getDenominator();
		if(den.compareTo(BigDecimal.ZERO) == 0) return toLog10(den);
		BigDecimal ans = getNumerator().setScale(64, RoundingMode.HALF_UP).divide(den, RoundingMode.HALF_UP);
		score = ans;
		return score;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		String seq = "";
		for(int strand = 0; strand < pfs.size()-1; ++strand) {
			seq += (pfs.get(strand).getSequence()+" | ");
		}
		int ind = seq.lastIndexOf("|");
		if(ind>=0) seq = new StringBuilder(seq).replace(ind, ind+1,"").toString();
		seq = seq.trim();
		sb.append("Seq: "+seq+", ");
		
		sb.append(String.format("log10(score): %12e, ", getScoreLog10()));
		for(int strand=0;strand<pfs.size();++strand) {
			BigDecimal z = pfs.get(strand)==null ? BigDecimal.ZERO : 
				pfs.get(strand).getZ();
			BigInteger numConfs = pfs.get(strand)==null ? BigInteger.ZERO :
				pfs.get(strand).getNumConfs();
			sb.append(String.format("pf: %2d, z: %12e, confs: %d ", strand, z, numConfs));
		}
		
		String ans = sb.toString().trim();
		return ans;//.substring(0,ans.length()-1);
	}
	
}
