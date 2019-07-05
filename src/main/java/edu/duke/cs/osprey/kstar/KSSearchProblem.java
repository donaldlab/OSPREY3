/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.kstar;

import java.math.BigInteger;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.emat.ReducedEnergyMatrix;
import edu.duke.cs.osprey.kstar.pruning.InvertedPruningMatrix;
import edu.duke.cs.osprey.kstar.pruning.ReducedPruningMatrix;
import edu.duke.cs.osprey.kstar.pruning.UnprunedPruningMatrix;
import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class KSSearchProblem extends SearchProblem {

	private static final long serialVersionUID = -6946555904198558059L;

	public ReducedPruningMatrix reducedMat = null;
	public InvertedPruningMatrix inverseMat = null;
	public ArrayList<String> flexibleRes = null;
	public ArrayList<ArrayList<String>> allowedAAs = null;
	public ArrayList<ArrayList<String>> reducedAllowedAAs = null;
	public ArrayList<Integer> posNums = null;
	public ParamSet params = null;


	public KSSearchProblem(ParamSet params, String name, String PDBFile, ArrayList<String> flexibleRes,
			ArrayList<ArrayList<String>> allowedAAs, boolean addWT, boolean contSCFlex, boolean useEPIC,
			EPICSettings epicSettings, boolean useTupExp, LUTESettings luteSettings,
			DEEPerSettings dset, ArrayList<String[]> moveableStrands,
			ArrayList<String[]> freeBBZones, boolean useEllipses, boolean useERef, boolean addResEntropy,
			boolean addWTRots, ResidueTermini termini, boolean useVoxelG) {

		super(name, PDBFile, flexibleRes, allowedAAs, addWT, contSCFlex, useEPIC, epicSettings, useTupExp, luteSettings,
				dset, moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWTRots, termini, useVoxelG, new ArrayList<>());

		this.params = params;
		
		this.allowedAAs = allowedAAs;
		this.reducedAllowedAAs = allowedAAs;
		this.posNums = getMaxPosNums();
	}


	public KSSearchProblem(SearchProblem sp1) {
		super(sp1);
	}


	public KSSearchProblem(KSSearchProblem other, 
			String newSPName, 
			ArrayList<ArrayList<String>> reducedAllowedAAs, 
			ArrayList<String> newFlexibleRes,
			ArrayList<Integer> newPosNums) {

		super(other);

		name = newSPName;
		this.allowedAAs = other.allowedAAs;
		this.reducedAllowedAAs = reducedAllowedAAs;		
		this.flexibleRes = newFlexibleRes;
		this.posNums = newPosNums;
		this.contSCFlex = other.contSCFlex;
		this.useERef = other.useERef;
		this.addResEntropy = other.addResEntropy;
		this.emat = other.emat;
		this.epicMat = other.epicMat;
		this.tupExpEMat = other.tupExpEMat;
		this.shellResidues = other.shellResidues;
		this.fullConfE = other.fullConfE;
		this.pruneMat = other.pruneMat;
		this.competitorPruneMat = other.competitorPruneMat;
		this.confSpace = other.confSpace;
	}

	public enum MatrixType {
		EMAT, EPICMAT;
	}

	public MatrixType getMatrixType() {
		//If the LUTE matrix is supposed to be used, it will be set as emat
		if(useEPIC) throw new UnsupportedOperationException("ERROR: EPIC is currently not supported in K*");

		return MatrixType.EMAT;
	}

	public EnergyMatrix getEnergyMatrix() {
		return emat;
		//If the LUTE matrix is supposed to be used, it will be set as emat
	}

	public String getMatrixFileName(MatrixType type) {
		return name + "." + type.name() + ".dat";
	}

	public String getEnergyMatrixFileName() {
		return name + ".EMAT.dat";
	}

	public ArrayList<Integer> getMaxPosNums() {
		ArrayList<Integer> ans = new ArrayList<>(allowedAAs.size());
		for(int i = 0; i < allowedAAs.size(); ++i) ans.add(i);
		return ans;
	}


	public InvertedPruningMatrix getInvertedFromUnreducedPruningMatrix(KSSearchProblem sp) {
		ReducedPruningMatrix unreduced = new ReducedPruningMatrix(this);

		// invert the QStar pruneMat
		InvertedPruningMatrix ans = new InvertedPruningMatrix(sp, unreduced.getUpdatedPruningMatrix());

		// throw exception if there is a sequence mismatch
		if(!pruneMatIsValid(ans))
			throw new RuntimeException("ERROR: pruning did not reduce RCs to sequence space of allowedAAs");

		BigInteger numConfs = numConfs(ans);
		return numConfs.compareTo(BigInteger.ZERO) == 0 ? null : ans;
	}


	public InvertedPruningMatrix getInvertedFromReducedPruningMatrix(KSSearchProblem sp) {

		// invert the QStar pruneMat
		InvertedPruningMatrix ans = new InvertedPruningMatrix(sp,
				((ReducedPruningMatrix)sp.reducedMat).getUpdatedPruningMatrix());

		// throw exception if there is a sequence mismatch
		if(!pruneMatIsValid(ans))
			throw new RuntimeException("ERROR: pruning did not reduce RCs to sequence space of allowedAAs");

		BigInteger numConfs = numConfs(ans);
		return numConfs.compareTo(BigInteger.ZERO) == 0 ? null : ans;
	}


	public UnprunedPruningMatrix getUnprunedPruningMatrix( KSSearchProblem sp, double pruningInterval ) {

		// invert the QStar pruneMat
		UnprunedPruningMatrix ans = new UnprunedPruningMatrix( sp, sp.reducedMat.getUpdatedPruningMatrix(), pruningInterval );

		// throw exception if there is a sequence mismatch
		if(!pruneMatIsValid(ans))
			throw new RuntimeException("ERROR: pruning did not reduce RCs to sequence space of allowedAAs");

		BigInteger numConfs = numConfs(ans);
		return numConfs.compareTo(BigInteger.ZERO) == 0 ? null : ans;
	}


	protected ReducedPruningMatrix getReducedPruningMatrix( SearchProblem sp ) {
		// see comets tree.dochildpruning

		// updated pruning matrix consists of rcs from this sequence only
		ReducedPruningMatrix ans = new ReducedPruningMatrix(this);

		// 1
		// prune all residues for other AA types at positions corresponding to this sequence
		for( int pos : posNums ) {
			for( int rc : sp.pruneMat.unprunedRCsAtPos(pos) ) {
				String rcAAType = confSpace.posFlex.get(pos).RCs.get(rc).AAType;
				if( !reducedAllowedAAs.get(posNums.indexOf(pos)).contains(rcAAType) ) 
					ans.getUpdatedPruningMatrix().markAsPruned(new RCTuple(pos,rc));
			}
		}

		// 2
		// prune all at positions not corresponding to this sequence
		for( int pos = 0; pos < sp.pruneMat.getNumPos(); ++pos ) {
			if(posNums.contains(pos)) continue;
			for( int rc : sp.pruneMat.unprunedRCsAtPos(pos) ) 
				ans.getUpdatedPruningMatrix().markAsPruned(new RCTuple(pos,rc));
		}

		// single sequence type dependent pruning for better efficiency
		//now do any consequent singles & pairs pruning
		int numUpdates = ans.countUpdates();
		int oldNumUpdates;

		double stericThresh = params == null ? 100.0 : params.getDouble("StericThresh");
		Pruner dee = new Pruner(sp, ans.getUpdatedPruningMatrix(), true, stericThresh, ans.getPruningInterval(), sp.useEPIC, sp.useTupExpForSearch);
		dee.setVerbose(false);

		do {//repeat as long as we're pruning things
			oldNumUpdates = numUpdates;
			dee.prune("GOLDSTEIN");
			dee.prune("GOLDSTEIN PAIRS FULL");
			numUpdates = ans.countUpdates();
		} while (numUpdates > oldNumUpdates && numConfs(ans).compareTo(BigInteger.valueOf(50000)) > 0);


		// throw exception if there is a sequence mismatch
		if(!pruneMatIsValid(ans))
			throw new RuntimeException("ERROR: pruning did not reduce RCs to sequence space of allowedAAs");

		BigInteger numConfs = numConfs(ans);
		return numConfs.compareTo(BigInteger.ZERO) == 0 ? null : ans;
	}


	public KSSearchProblem getReducedSearchProblem( String name, 
			ArrayList<ArrayList<String>> allowedAAs, 
			ArrayList<String> flexRes, 
			ArrayList<Integer> posNums ) {

		// Create a version of the search problem restricted to the specified sequence (list of amino acid names)
		// the constructor creates a new confspace object
		KSSearchProblem reducedSP = new KSSearchProblem(this, name, allowedAAs, flexRes, posNums);

		reducedSP.reducedMat = reducedSP.getReducedPruningMatrix(reducedSP); // for q*, q'

		if(reducedSP.reducedMat != null)
			reducedSP.inverseMat = reducedSP.getInvertedFromReducedPruningMatrix(reducedSP); // for p*

		return reducedSP;
	}


	public ArrayList<String> getAAsAtPos( PruningMatrix pruneMat, int pos ) {
		ArrayList<String> ans = new ArrayList<>();

		ArrayList<Integer> rcsAtPos = pruneMat.unprunedRCsAtPos(pos);
		for( int RCNum : rcsAtPos ) {
			int pos1 = posNums.get(pos);
			String AAType = confSpace.posFlex.get(pos1).RCs.get(RCNum).AAType;
			if(!ans.contains(AAType)) ans.add(AAType);
		}

		return ans;
	}


	public ArrayList<Integer> rcsAtPos( PruningMatrix pruneMat, int pos, String aaType, boolean pruned ) {
		ArrayList<Integer> ans = new ArrayList<>();

		ArrayList<Integer> rcsAtPos = pruned ? pruneMat.prunedRCsAtPos(pos) : pruneMat.unprunedRCsAtPos(pos);
		for( int RCNum : rcsAtPos ) {
			int pos1 = posNums.get(pos);
			String AAType = confSpace.posFlex.get(pos1).RCs.get(RCNum).AAType;
			if(!AAType.equalsIgnoreCase(aaType)) continue;
			ans.add(RCNum);
		}

		return ans;
	}


	private ArrayList<ArrayList<String>> getAAsAtPos( PruningMatrix pruneMat ) {
		ArrayList<ArrayList<String>> ans = new ArrayList<>();

		for( int pos = 0; pos < pruneMat.getNumPos(); ++pos )
			ans.add(getAAsAtPos(pruneMat, pos));

		return ans;
	}


	private boolean pruneMatIsValid( PruningMatrix pruneMat ) {

		ArrayList<ArrayList<String>> pruneMatAAs = getAAsAtPos(pruneMat);

		if(pruneMatAAs.size() != reducedAllowedAAs.size()) 
			return false;

		for(int pos = 0; pos < reducedAllowedAAs.size(); ++pos) {
			for(String aaAtPos : pruneMatAAs.get(pos)) {
				if(!reducedAllowedAAs.get(pos).contains(aaAtPos))
					return false;
			}
		}

		return true;
	}


	public BigInteger numConfs( PruningMatrix pruneMat ) {

		if( pruneMat == null ) return BigInteger.ZERO;

		BigInteger ans = BigInteger.ONE;

		for( int pos = 0; pos < pruneMat.getNumPos(); ++pos ) {
			long numRCs = pruneMat.unprunedRCsAtPos(pos).size();
			if(numRCs == 0) return BigInteger.ZERO;
			ans = ans.multiply( BigInteger.valueOf( numRCs ) );
		}

		return ans;
	}


	public EnergyMatrix getReducedEnergyMatrix() {
		// i only use this for partial sequences
		if(posNums.size() == confSpace.numPos)
			return emat;

		return new ReducedEnergyMatrix(this, emat);
	}


	public MultiTermEnergyFunction decompMinimizedEnergy(int[] conf){
		//Minimized energy of the conformation
		//whose RCs are listed for all flexible positions in conf
		MultiTermEnergyFunction mef = confSpace.getDecomposedMinimizedEnergy(conf, fullConfE, null);

		double E = mef.getPreCompE();

		E += emat.getConstTerm();

		if(useERef)
			E -= emat.geteRefMat().confERef(conf);

		if(addResEntropy)
			E += confSpace.getConfResEntropy(conf);            

		mef.setPreCompE(E);

		return mef;
	}


	public void mergeResiduePositions(int... posToCombine) {

		EnergyMatrixCalculator emc = new EnergyMatrixCalculator(confSpace, 
				shellResidues, useEPIC, reducedMat, epicSettings, false, emat);

		emc.addEnergyTerms(false, posToCombine);
	}


	public double lowerBoundContribByRC(int pos, int[] conf, int numResInHot) {
		double bound = emat.rcContribAtPos(pos, conf, numResInHot);	
		return bound;
	}
}
