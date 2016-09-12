package edu.duke.cs.osprey.kstar;

import java.math.BigInteger;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.emat.ReducedEnergyMatrix;
import edu.duke.cs.osprey.kstar.pruning.InvertedPruningMatrix;
import edu.duke.cs.osprey.kstar.pruning.ReducedPruningMatrix;
import edu.duke.cs.osprey.kstar.pruning.UnprunedPruningMatrix;
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
	
	
	public KSSearchProblem(String name, String PDBFile, ArrayList<String> flexibleRes,
			ArrayList<ArrayList<String>> allowedAAs, boolean addWT, boolean contSCFlex, boolean useEPIC,
			EPICSettings epicSettings, boolean useTupExp, LUTESettings luteSettings,
            DEEPerSettings dset, ArrayList<String[]> moveableStrands,
			ArrayList<String[]> freeBBZones, boolean useEllipses, boolean useERef, boolean addResEntropy,
			boolean addWTRots, KSTermini termini, boolean useVoxelG) {
		
		super(name, PDBFile, flexibleRes, allowedAAs, addWT, contSCFlex, useEPIC, epicSettings, useTupExp, luteSettings,
                        dset, moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWTRots, termini, useVoxelG);
		
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
	
	
	public void loadMatrix() {
		MatrixType type = getMatrixType();
		
		switch (type) {

		case EMAT: 
			loadEnergyMatrix();
			break;

		case TUPEXPEMAT: 
			loadTupExpEMatrix();
			break;

		case EPICMAT:
			loadEPICMatrix();
			break;

		default:	
			throw new RuntimeException("ERROR: unsupported energy matrix type: " + type);
		}
	}
	
	
	public MatrixType getMatrixType() {
		if(useTupExpForSearch) return MatrixType.TUPEXPEMAT;

		else if(useEPIC) return MatrixType.EPICMAT;

		return MatrixType.EMAT;
	}


	public EnergyMatrix getEnergyMatrix() {

		MatrixType type = getMatrixType();
		switch (type) {

		case EMAT: 
			return emat;

		case TUPEXPEMAT: 
			return tupExpEMat;

		default:
			throw new RuntimeException("ERROR: AAO has not added support for type " + type);
		}
	}

	
	public void setEnergyMatrix(EnergyMatrix e) {

		MatrixType type = getMatrixType();
		switch (type) {

		case EMAT: 
			emat = e;
			break;

		case TUPEXPEMAT: 
			tupExpEMat = e;
			break;

		default:
			throw new RuntimeException("ERROR: AAO has not added support for type " + type);
		}
	}

	
	
	public String getMatrixFileName(MatrixType type) {
		return name + "." + type.name() + ".dat";
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
