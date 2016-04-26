package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.structure.Residue;

public class RCEnergyContribs {

	private PFAbstract pf = null;
	private HashMap<Integer, Double> minEContribs = null;
	private HashMap<Integer, Double> pwEContribs = null;
	private double minE = Double.POSITIVE_INFINITY;
	private double pwE = Double.POSITIVE_INFINITY;
	private HashSet<String> flexRes = null;
	ArrayList<BoundErrorByPos> boundErrorByPos = null;
	
	
	public RCEnergyContribs( PFAbstract pf, MultiTermEnergyFunction mef, int[] conf ) {
		this.pf = pf;
		
		flexRes = new HashSet<>();
		flexRes.addAll(AllowedSeqs.getFlexResFromSeq(pf.getSequence()));
		
		computeEContrib(mef, conf);
		
		boundErrorByPos = getBoundErrorByPos();
	}

	
	protected class BoundErrorByPos {
		
		int pos;
		double boundError;
		
		public BoundErrorByPos(int pos, double boundError) {
			this.pos = pos;
			this.boundError = boundError;
		}
	}
	
	
	protected class BoundErrorByPosComparator implements Comparator<BoundErrorByPos> {
		@Override
		public int compare(BoundErrorByPos o1, BoundErrorByPos o2) {
			if(o1.boundError < o2.boundError) return -1;
			return 1;
		}
		
	}
	
	
	protected ArrayList<BoundErrorByPos> getBoundErrorByPos() {
		
		if( boundErrorByPos != null ) return boundErrorByPos;
		
		int numPos = pf.getSearchProblem().confSpace.numPos;
		
		ArrayList<BoundErrorByPos> boundErrors = new ArrayList<>(numPos);
		for( int i = 0; i < numPos; ++i ) {
			//boundErrors.add( new BoundErrorByPos(i, Math.abs(pwEContribs.get(i)-minEContribs.get(i))) );
			boundErrors.add( new BoundErrorByPos(i, pwEContribs.get(i)-minEContribs.get(i)) );
		}
		
		// filter positions that already participate in tripples
		for( Iterator<BoundErrorByPos> iterator = boundErrors.iterator(); iterator.hasNext(); ) {
			BoundErrorByPos errorAtPos = iterator.next();
			if(pf.getPosInTripple().contains(errorAtPos.pos))
				iterator.remove();
		}
		
		// sort bound bound errors by contribution
		Collections.sort(boundErrors, new BoundErrorByPosComparator());
		return boundErrors;
	}
	
	
	// get bound error attributed to the top num rcs, sorted by error contribution
	public double getPercentBoundError( int num ) {
		
		if(num < 3) throw new RuntimeException("ERROR: num must be >= 3");
		
		if( num > boundErrorByPos.size() ) return 0;
		
		double errorByNum = 0;
		for(int i = 0; i < num; ++i) errorByNum += boundErrorByPos.get(i).boundError;
		
		double totalBoundError = getBoundError();
		errorByNum /= totalBoundError;
		
		return errorByNum;
	}
	
	
	public ArrayList<Integer> getPosCausingError( int num ) {
		
		if(num < 3) throw new RuntimeException("ERROR: num must be >= 3");
		
		if(num > boundErrorByPos.size()) return null;
		
		ArrayList<Integer> ans = new ArrayList<>(num);
		for(int i = 0; i < num; ++i) ans.add(boundErrorByPos.get(i).pos);
		return ans;
	}
	
	
	public double getBoundError() {
		//return Math.abs(pwE-minE);
		return pwE-minE;
	}
	
	
	protected void computeEContrib( MultiTermEnergyFunction mef, int[] conf ) {
		
		minE = getMinE(mef, true);
		pwE = getPWLB(conf, true);
		
		if(minE < pwE) throw new RuntimeException("ERROR: minE: " + minE + " must be >= pairwiseE: " + pwE);
	}
	
	
	public double getMinE( MultiTermEnergyFunction mef, boolean addTemplateSelfE ) {
		
		if(minE != Double.POSITIVE_INFINITY) return minE;
		
		minEContribs = getMinEContrib(mef);
		minE = getE(minEContribs);
		if(addTemplateSelfE) minE += getConstTerm();
		return minE;
	}
	
	
	public double getPWLB( int[] conf, boolean addTemplateSelfE ) {
		
		if(pwE != Double.POSITIVE_INFINITY) return pwE;
		
		pwEContribs = getPairwiseEContrib(conf);
		pwE = getE(pwEContribs);
		if(addTemplateSelfE) pwE += getConstTerm();
		return pwE;
	}
	
	
	private double getConstTerm() {
		return pf.getSearchProblem().getConstTerm();
	}
	
	
	private double getE( HashMap<Integer, Double> map ) {
		double E = 0;
		for( double partialE : map.values() )  E += partialE;
		return E;
	}
	

	// energy contributions exclude reference energy and template self-energy
	private HashMap<Integer, Double> getMinEContrib( MultiTermEnergyFunction mef ) {
		
		int numPos = pf.getSearchProblem().confSpace.numPos;
		HashMap<Integer, Double> ans = new HashMap<>(numPos);

		for( int i = 0; i < numPos; ++i ) {
			double intra = 0;
			double pairwise = 0;

			// get intra contribution
			String aa = pf.getSequence().get(i).split("-")[0].toLowerCase();
			String resNum = pf.getSequence().get(i).split("-")[1].toLowerCase();

			for( int j = 0; j < mef.getTerms().size(); ++j ) {

				if( mef.getTerms().get(j) instanceof SingleResEnergy ) {
					Residue res = ((SingleResEnergy)mef.getTerms().get(j)).getRes();

					if(res.fullName.toLowerCase().contains(aa) && res.fullName.contains(resNum))
						intra += mef.getCoeffs().get(j) * mef.getTerms().get(j).getEnergy();
				}

				else if( mef.getTerms().get(j) instanceof ResPairEnergy ) {
					Residue res1 = ((ResPairEnergy)mef.getTerms().get(j)).getRes1();
					Residue res2 = ((ResPairEnergy)mef.getTerms().get(j)).getRes2();

					if( (res1.fullName.toLowerCase().contains(aa) && res1.fullName.contains(resNum)) 
							|| (res2.fullName.toLowerCase().contains(aa) && res2.fullName.contains(resNum)) ) {
					
						double partialPW = mef.getCoeffs().get(j) * mef.getTerms().get(j).getEnergy();
						
						if( flexRes.contains(res1.getPDBResNumber()) && flexRes.contains(res2.getPDBResNumber()) )
							partialPW *= 0.5;
						
						pairwise += partialPW;
					
					}
				}
			}

			ans.put(i, intra + pairwise);
		}

		return ans;
	}


	private HashMap<Integer, Double> getPairwiseEContrib( int[] conf ) {
		
		int numPos = pf.getSearchProblem().confSpace.numPos;
		HashMap<Integer, Double> ans = new HashMap<>(numPos);

		for( int pos = 0; pos < conf.length; ++pos ) {
			double E = pf.getSearchProblem().lowerBoundContribByRC(pos, conf);
			ans.put(pos, E);
		}

		return ans;
	}

}
