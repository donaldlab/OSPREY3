package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;

import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.structure.Residue;

public class RCEnergyContribs {

	private PFAbstract pf = null;
	private HashMap<Integer, Double> flexResMinEContribs = null;
	private HashMap<Integer, Double> flexResPWLBContribs = null;
	private double flexResMinE = Double.POSITIVE_INFINITY;
	private double flexResPWLB = Double.POSITIVE_INFINITY;
	private ArrayList<BoundErrorByPos> boundErrorByPos = null;
	private HashMap<String, String> resNum2AA = null;
	private double minE = 0.0;
	private double pwLB = 0.0;


	public RCEnergyContribs( PFAbstract pf, MultiTermEnergyFunction mef, int[] conf ) {
		this.pf = pf;

		resNum2AA = makeResNum2AA();

		minE = mef.getPreCompE();
		pwLB = pf.getConfBound(null, conf);

		computeEContrib(mef, conf);

		boundErrorByPos = getBoundErrorByPos();
	}


	private HashMap<String, String> makeResNum2AA() {

		HashMap<String, String> ans = new HashMap<>(pf.getSequence().size());

		for( String res : pf.getSequence() ) {
			// get intra contribution
			String aa = res.split("-")[0];
			String resNum = res.split("-")[1];

			ans.put(resNum, aa);
		}

		return ans;
	}

	
	protected boolean isRes( Residue res, String resNum, String aa ) {
		return res.getPDBResNumber().compareTo(resNum) == 0 && res.fullName.contains(aa);
	}

	
	protected boolean inFlexRes( Residue res ) {

		String aa = resNum2AA.get(res.getPDBResNumber());
		if(aa == null) return false;

		if(!res.fullName.contains(aa)) return false;

		return true;
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

		int numPos = pf.getReducedSearchProblem().confSpace.numPos;

		ArrayList<BoundErrorByPos> boundErrors = new ArrayList<>(numPos);
		for( int pos = 0; pos < numPos; ++pos ) {
			//boundErrors.add( new BoundErrorByPos(i, Math.abs(pwEContribs.get(i)-minEContribs.get(i))) );
			boundErrors.add( new BoundErrorByPos(pos, flexResPWLBContribs.get(pos)-flexResMinEContribs.get(pos)) );
		}

		// filter positions that already participate in tripples
		for( Iterator<BoundErrorByPos> iterator = boundErrors.iterator(); iterator.hasNext(); ) {
			BoundErrorByPos errorAtPos = iterator.next();
			if(pf.HOTsContains(errorAtPos.pos))
				iterator.remove();
		}

		// sort bound bound errors by contribution
		Collections.sort(boundErrors, new BoundErrorByPosComparator());
		return boundErrors;
	}


	// get bound error attributed to the top num rcs, sorted by error contribution
	public double getPercentErrorForTopPos( int num ) {

		if(num < 3) throw new RuntimeException("ERROR: num must be >= 3");

		if( num > boundErrorByPos.size() ) return 0;

		double percentErrorByNum = 0;
		for(int pos = 0; pos < num; ++pos) percentErrorByNum += boundErrorByPos.get(pos).boundError;

		double allFlexResBoundError = getReComputedBoundError();
		percentErrorByNum /= allFlexResBoundError;

		//if(percentErrorByNum > 1.0) 
		//	throw new RuntimeException("ERROR: percent error due to top " + num + " rotamers: " + percentErrorByNum + " is > 1");

		return percentErrorByNum;
	}


	public int[] getTopPosCausingError( int num ) {

		if(num < 3) throw new RuntimeException("ERROR: num must be >= 3");

		int[] ans = new int[num];

		if(num > boundErrorByPos.size()) return ans;

		for(int i = 0; i < num; ++i) ans[i] = boundErrorByPos.get(i).pos;

		Arrays.sort(ans);

		return ans;
	}


	public double getReComputedBoundError() {
		return flexResPWLB-flexResMinE;
	}
	
	
	public double getPercentBoundError() {
		return getReComputedBoundError()/flexResMinE;
	}


	public double getBoundError() {
		return pwLB - minE;
	}


	protected void computeEContrib( MultiTermEnergyFunction mef, int[] conf ) {

		flexResMinE = getFlexResMinE(mef, true);
		flexResPWLB = getFlexResPWLB(conf, true);

		// error checks
		if(flexResMinE < flexResPWLB) throw new RuntimeException("ERROR: flexResMinE: " + flexResMinE + " must be >= flexResPairWiseLB: " + flexResPWLB);
		if(Math.abs(flexResPWLB-pwLB) > 0.0001) throw new RuntimeException("ERROR: re-computed pwLB: " + flexResPWLB + " != actual pwLB: " + pwLB);
		if(Math.abs(flexResMinE-minE) > 0.0001) throw new RuntimeException("ERROR: re-computed minE: " + flexResPWLB + " != actual minE: " + minE);
	}


	public double getFlexResMinE( MultiTermEnergyFunction mef, boolean addTemplateSelfE ) {

		if(flexResMinE != Double.POSITIVE_INFINITY) return flexResMinE;

		flexResMinEContribs = getFlexResMinEContribs(mef);
		flexResMinE = sumETerms(flexResMinEContribs);
		if(addTemplateSelfE) flexResMinE += getConstTerm();
		return flexResMinE;
	}


	public double getFlexResPWLB( int[] conf, boolean addTemplateSelfE ) {

		if(flexResPWLB != Double.POSITIVE_INFINITY) return flexResPWLB;

		flexResPWLBContribs = getFlexResPWLBContribs(conf);
		flexResPWLB = sumETerms(flexResPWLBContribs);
		if(addTemplateSelfE) flexResPWLB += getConstTerm();
		return flexResPWLB;
	}


	private double getConstTerm() {
		return pf.getReducedSearchProblem().getEnergyMatrix().getConstTerm();
	}


	private double sumETerms( HashMap<Integer, Double> map ) {
		double E = 0;
		for( double partialE : map.values() )  E += partialE;
		return E;
	}


	// store single and pairwise energy contributions involving flexible residues only.
	private HashMap<Integer, Double> getFlexResMinEContribs( MultiTermEnergyFunction mef ) {

		int numPos = pf.getReducedSearchProblem().confSpace.numPos;
		HashMap<Integer, Double> ans = new HashMap<>(numPos);

		for( int pos = 0; pos < numPos; ++pos ) {
			double intra = 0;
			double pairwise = 0;

			// get intra contribution
			String aa = pf.getSequence().get(pos).split("-")[0];
			String resNum = pf.getSequence().get(pos).split("-")[1];

			for( int j = 0; j < mef.getTerms().size(); ++j ) {
				
				if( mef.getTerms().get(j) instanceof SingleResEnergy ) {
					Residue res = ((SingleResEnergy)mef.getTerms().get(j)).getRes();
					
					if(isRes(res, resNum, aa))
						intra += mef.getCoeffs().get(j) * mef.getTerms().get(j).getEnergy();
				}

				else if( mef.getTerms().get(j) instanceof ResPairEnergy ) {
					Residue res1 = ((ResPairEnergy)mef.getTerms().get(j)).getRes1();
					Residue res2 = ((ResPairEnergy)mef.getTerms().get(j)).getRes2();
					
					if( (isRes(res1, resNum, aa) && inFlexRes(res2)) || (isRes(res2, resNum, aa) && inFlexRes(res1)) ) {
						double partialPW = mef.getCoeffs().get(j) * mef.getTerms().get(j).getEnergy();
						pairwise += 0.5 * partialPW;
					}
					
					// See EnergyFunctionGenerator.intraAndShellEnergy(...)
					// intra energy consists of single residue energy and sum of pairs consisting of residue to shell energies
					else if( (isRes(res1, resNum, aa) && !inFlexRes(res2)) || (isRes(res2, resNum, aa) && !inFlexRes(res1)) ) {						
						intra += mef.getCoeffs().get(j) * mef.getTerms().get(j).getEnergy();
					}
				}
			}

			ans.put(pos, intra + pairwise);
		}

		return ans;
	}


	private HashMap<Integer, Double> getFlexResPWLBContribs( int[] conf ) {

		int numPos = pf.getReducedSearchProblem().confSpace.numPos;
		HashMap<Integer, Double> ans = new HashMap<>(numPos);

		for( int pos = 0; pos < conf.length; ++pos ) {
			double E = pf.getReducedSearchProblem().lowerBoundContribByRC(pos, conf, PFAbstract.getHotNumRes());
			ans.put(pos, E);
		}

		return ans;
	}

}
