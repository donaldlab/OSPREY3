/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.Termini;
import edu.duke.cs.osprey.kstar.emat.ReducedEnergyMatrix;
import edu.duke.cs.osprey.kstar.pruning.InvertedPruningMatrix;
import edu.duke.cs.osprey.kstar.pruning.ReducedPruningMatrix;
import edu.duke.cs.osprey.kstar.pruning.UnprunedPruningMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.ConfETupleExpander;
import edu.duke.cs.osprey.tupexp.TupExpChooser;
import java.io.Serializable;
import java.math.BigInteger;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 * @author Adegoke Ojewole (ao68@duke.edu)
 */
@SuppressWarnings("serial")
public class SearchProblem implements Serializable {
	//This object keeps track of the positions and the possible assignments for them, as used in the search algorithms
	//generally these will be rDesidues (or super-residues) and their RCs; subclass SearchProblem to change this

	//We keep a ConfSpace together with "annotations" that help us find its GMEC, partition functions, etc.
	//annotations are based on RCs and indicate pairwise energies, pruning information, etc.
	//they also include iterators over RCs and pairs of interest

	public ConfSpace confSpace;

	public EnergyMatrix emat = null;//energy matrix.  Meanings:
	//-Defines full energy in the rigid case
	//-Defines lower bound in the continuous case
	//-emat + epicm = full energy if using EPIC for search 

	public EPICMatrix epicMat = null;//EPIC matrix, to be used if appropriate
	public EPICSettings epicSettings = null;

	public EnergyMatrix tupExpEMat;//Defines full energy in the continuous, tuple-expander case

	public EnergyFunction fullConfE;//full energy for any conformation
	public ArrayList<Residue> shellResidues;//non-flexible residues to be accounted for in energy calculations

	public String name;//a human-readable name, which will also be used to name stored energy matrices, etc.

	public PruningMatrix pruneMat = null;
	public ReducedPruningMatrix reducedMat = null;
	public InvertedPruningMatrix inverseMat = null;

	public boolean contSCFlex;

	public ArrayList<String> flexibleRes = null;
	public ArrayList<ArrayList<String>> allowedAAs = null;
	public ArrayList<ArrayList<String>> reducedAllowedAAs = null;
	public ArrayList<Integer> posNums = null;

	public DEEPerSettings dset;
	public ArrayList<String[]> moveableStrands;
	public ArrayList<String[]> freeBBZones;

	public PruningMatrix competitorPruneMat;//a pruning matrix performed at pruning interval 0,
	//to decide which RC tuples are valid competitors for pruning

	public boolean useEPIC = false;
	public boolean useTupExpForSearch = false;//use a tuple expansion to approximate the energy as we search
	public boolean useEllipses = false;

	public boolean useERef = false;
	public boolean addResEntropy = false;
	public Termini limits = null;


	public SearchProblem(SearchProblem sp1){//shallow copy
		confSpace = sp1.confSpace;
		emat = sp1.emat;
		epicMat = sp1.epicMat;
		epicSettings = sp1.epicSettings;
		tupExpEMat = sp1.tupExpEMat;

		fullConfE = sp1.fullConfE;
		shellResidues = sp1.shellResidues;
		name = sp1.name + System.currentTimeMillis();//probably will want to change this to something more meaningful

		pruneMat = sp1.pruneMat;
		competitorPruneMat = sp1.competitorPruneMat;

		contSCFlex = sp1.contSCFlex;
		useEPIC = sp1.useEPIC;
		useTupExpForSearch = sp1.useTupExpForSearch;

		useERef = sp1.useERef;
		addResEntropy = sp1.addResEntropy;
		limits = sp1.limits;
	}


	public SearchProblem(SearchProblem origSP, 
			String newSPName, 
			ArrayList<ArrayList<String>> reducedAllowedAAs, 
			ArrayList<String> newFlexibleRes,
			ArrayList<Integer> newPosNums) {

		name = newSPName;

		// flexible residues are the same
		this.allowedAAs = origSP.allowedAAs;
		this.reducedAllowedAAs = reducedAllowedAAs;		
		flexibleRes = newFlexibleRes;
		this.posNums = newPosNums;

		contSCFlex = origSP.contSCFlex;
		useTupExpForSearch = origSP.useTupExpForSearch;
		useEllipses = origSP.useEllipses;
		useEPIC = origSP.useEPIC;
		epicSettings = origSP.epicSettings;

		useERef = origSP.useERef;
		addResEntropy = origSP.addResEntropy;

		emat = origSP.emat;
		epicMat = origSP.epicMat;
		tupExpEMat = origSP.tupExpEMat;

		dset = origSP.dset;
		moveableStrands = origSP.moveableStrands; 
		freeBBZones = origSP.freeBBZones;

		limits = origSP.limits;

		shellResidues = origSP.shellResidues;
		fullConfE = origSP.fullConfE;

		pruneMat = origSP.pruneMat;
		competitorPruneMat = origSP.competitorPruneMat;

		confSpace = origSP.confSpace;
	}


	public SearchProblem(String name, String PDBFile, ArrayList<String> flexibleRes, ArrayList<ArrayList<String>> allowedAAs, boolean addWT,
			boolean contSCFlex, boolean useEPIC, EPICSettings epicSettings, boolean useTupExp, DEEPerSettings dset, 
			ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, boolean useEllipses, boolean useERef,
			boolean addResEntropy, Termini limits){

		confSpace = new ConfSpace(PDBFile, flexibleRes, allowedAAs, addWT, contSCFlex, dset, moveableStrands, freeBBZones, useEllipses, limits);
		this.name = name;

		this.flexibleRes = flexibleRes;
		this.allowedAAs = allowedAAs;

		this.contSCFlex = contSCFlex;
		this.useTupExpForSearch = useTupExp;
		this.useEllipses = useEllipses;
		this.useEPIC = useEPIC;
		this.epicSettings = epicSettings;

		this.useERef = useERef;
		this.addResEntropy = addResEntropy;

		this.dset = dset;
		this.moveableStrands = moveableStrands;
		this.freeBBZones = freeBBZones;
		this.limits = limits;

		//energy function setup
		EnergyFunctionGenerator eGen = EnvironmentVars.curEFcnGenerator;
		decideShellResidues(eGen.distCutoff);
		fullConfE = eGen.fullConfEnergy(confSpace,shellResidues);

		this.reducedAllowedAAs = allowedAAs;
		this.posNums = getMaxPosNums();
	}


	public ArrayList<Integer> getMaxPosNums() {
		ArrayList<Integer> ans = new ArrayList<>(allowedAAs.size());
		for(int i = 0; i < allowedAAs.size(); ++i) ans.add(i);
		return ans;
	}

	public void mergeResiduePositions(int... posToCombine) {

		EnergyMatrixCalculator emc = new EnergyMatrixCalculator(confSpace, 
				shellResidues, useEPIC, pruneMat, epicSettings, false, emat);

		emc.addEnergyTerms(false, posToCombine);
	}


	private void decideShellResidues(double distCutoff){
		//Decide what non-flexible residues need to be accounted for in energy calculations

		ArrayList<Residue> flexibleResidues = new ArrayList<>();//array list for these so they stay in order
		for(PositionConfSpace pcs : confSpace.posFlex)
			flexibleResidues.add(pcs.res);

		//we'll decide what shell residues to include by using a simple distance cutoff with
		//the current conformation,
		//rather than doing a conformational search for the minimal distance (with respect to conformations)
		//of a shell residue to any flexible residues
		//the distance cutoff can be increased to accommodate this if desired.
		shellResidues = new ArrayList<>();

		for(Residue nonFlexRes : confSpace.m.residues){
			if(!flexibleResidues.contains(nonFlexRes)){//residue is not flexible

				for(Residue flexRes : flexibleResidues){
					double dist = flexRes.distanceTo(nonFlexRes);
					if(dist<=distCutoff){
						shellResidues.add(nonFlexRes);//close enough to flexRes that we should add it
						break;
					}
				}
			}
		}
	}


	public double minimizedEnergy(int[] conf){
		//Minimized energy of the conformation
		//whose RCs are listed for all flexible positions in conf
		double E = confSpace.minimizeEnergy(conf, fullConfE, null);

		E += getConstTerm();

		if(useERef)
			E -= emat.geteRefMat().confERef(conf);

		if(addResEntropy)
			E += confSpace.getConfResEntropy(conf);            

		return E;
	}


	public MultiTermEnergyFunction decompMinimizedEnergy(int[] conf){
		//Minimized energy of the conformation
		//whose RCs are listed for all flexible positions in conf
		MultiTermEnergyFunction mef = confSpace.getDecomposedMinimizedEnergy(conf, fullConfE, null);

		double E = mef.getPreCompE();

		E += getConstTerm();

		if(useERef)
			E -= emat.geteRefMat().confERef(conf);

		if(addResEntropy)
			E += confSpace.getConfResEntropy(conf);            

		mef.setPreCompE(E);

		return mef;
	}


	public void outputMinimizedStruct(int[] conf, String PDBFileName){
		//Output minimized conformation to specified file
		//RCs are listed for all flexible positions in conf
		//Note: eRef not included (just minimizing w/i voxel)
		confSpace.minimizeEnergy(conf, fullConfE, PDBFileName);
	}


	public double approxMinimizedEnergy(int[] conf){
		//EPIC or other approximation for the minimized energy of the conformation
		//whose RCs are listed for all flexible positions in conf

		if( useTupExpForSearch ){//use tup-exp E-matrix direectly
			return tupExpEMat.confE(conf);
		}
		else if( useEPIC ){//EPIC w/o tup-exp
			return EPICMinimizedEnergy(conf);
		}

		else
			throw new RuntimeException("ERROR: Asking searchSpace to approximate minimized energy but using a non-approximating method");
	}


	public double EPICMinimizedEnergy(int[] conf){
		//approximate energy using EPIC
		double bound = emat.confE(conf);//emat contains the pairwise lower bounds
		double contPart = epicMat.minContE(conf);
		//EPIC handles the continuous part (energy - pairwise lower bounds)

		return bound+contPart;
	}



	public double lowerBound(int[] conf){
		//Argument: RC assignments for all the flexible residues (RCs defined in resFlex)
		//return lower bound on energy for the conformational space defined by these assignments
		//based on precomputed energy matrix (including EPIC if indicated)

		double bound = emat.confE(conf);//the energy recorded by the matrix is 
		//the pairwise lower bounds
		return bound;
	}


	public double lowerBoundContribByRC(int pos, int[] conf, int numResInHot) {
		double bound = emat.rcContribAtPos(pos, conf, numResInHot);	
		return bound;
	}


	public double getConstTerm() {
		return emat.getConstTerm();
	}


	//LOADING AND PRECOMPUTATION OF ENERGY MATRIX-TYPE OBJECTS (regular energy matrix, tup-exp and EPIC matrices)

	// AAO
	public MatrixType getMatrixType() {
		if(useTupExpForSearch) return MatrixType.TUPEXPEMAT;

		// skipping epic for now
		// else if(useEPIC) return MatrixType.EPICMAT;

		return MatrixType.EMAT;
	}

	// AAO
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

	// AAO
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

	// AAO it's easier to specify the matrix type
	public void loadEnergyMatrix(MatrixType type) {
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

	public void loadEnergyMatrix(){
		loadMatrix(MatrixType.EMAT);
	}

	public void loadTupExpEMatrix(){
		loadMatrix(MatrixType.TUPEXPEMAT);
	}

	public void loadEPICMatrix(){
		loadMatrix(MatrixType.EPICMAT);
	}

	public enum MatrixType {
		EMAT, TUPEXPEMAT, EPICMAT;
	}


	public String getMatrixFileName(MatrixType type) {
		return name + "." + type.name() + ".dat";
	}


	//load the specified matrix; if the right file isn't available then compute and store it
	private void loadMatrix(MatrixType type){

		String matrixFileName = getMatrixFileName(type);
		//matrix file names are determined by the name of the search problem

		if(!loadMatrixFromFile( type, matrixFileName )){
			TupleMatrix matrix = calcMatrix(type);
			ObjectIO.writeObject( matrix, matrixFileName );
			loadMatrixFromFile( type, matrixFileName );
		}
	}


	//compute the matrix of the specified type
	private TupleMatrix calcMatrix(MatrixType type){

		if(type == MatrixType.EMAT){
			EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpace,shellResidues,
					useERef,addResEntropy);

			emCalc.calcPEM();
			return emCalc.getEMatrix();
		}
		else if(type == MatrixType.EPICMAT){
			EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpace,shellResidues,
					pruneMat,epicSettings);
			emCalc.calcPEM();
			return emCalc.getEPICMatrix();
		}
		else {
			//need to calculate a tuple-expansion matrix

			double errorThresh = 0.01;

			ConfETupleExpander expander = new ConfETupleExpander(this);//make a tuple expander
			TupleEnumerator tupEnum = new TupleEnumerator(pruneMat,emat,confSpace.numPos);
			TupExpChooser chooser = new TupExpChooser(expander, tupEnum);//make a chooser to choose what tuples will be in the expansion

			double curResid = chooser.calcPairwiseExpansion();//start simple...

			if(curResid > errorThresh){//go to triples if needed
				System.out.println("EXPANDING PAIRWISE EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (2 PARTNERS)...");
				curResid = chooser.calcExpansionResTriples(2);
			}
			if(curResid > errorThresh){//go to 5 partners if still need better resid...
				System.out.println("EXPANDING EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (5 PARTNERS)...");
				curResid = chooser.calcExpansionResTriples(5);
			}
			if(curResid > errorThresh){
				System.out.println("WARNING: Desired LUTE residual threshold "+
						errorThresh+" not reached; best="+curResid);
			}

			return expander.getEnergyMatrix();//get the final energy matrix from the chosen expansion
		}
	}



	boolean loadMatrixFromFile(MatrixType type, String matrixFileName){
		//try loading the specified matrix from a file
		//return true if successful, false if not, in which case we'll have to compute it
		//also if the matrix's pruning interval is too low, it may be missing some RCs
		//that are unpruned at our current pruningInterval, so we have to recompute
		Object matrixFromFile = ObjectIO.readObject(matrixFileName, true);

		if(type == MatrixType.EMAT)
			emat = (EnergyMatrix) matrixFromFile;
		else if(type == MatrixType.EPICMAT)
			epicMat = (EPICMatrix) matrixFromFile;
		else //tup-exp
			tupExpEMat = (EnergyMatrix) matrixFromFile;

		if(matrixFromFile==null)//unsuccessful loading leaves null emat
			return false;


		//check pruning interval.  Current interval is in pruneMat if we have pruned already;
		//if not then we need a matrix with infinite pruning interval (valid for all RCs).
		double matrixPruningInterval = ((TupleMatrix)matrixFromFile).getPruningInterval();

		if( matrixPruningInterval == Double.POSITIVE_INFINITY )//definitely valid
			return true;
		else {
			//excludes some RCs...check against pruneMat pruning interval
			if(pruneMat==null){
				throw new RuntimeException("ERROR: Trying to load pruning-dependent tuple matrix"
						+ "(EPIC or tup-exp) but haven't pruned yet");
			}

			return ( matrixPruningInterval >= pruneMat.getPruningInterval() );
		}
	}


	public InvertedPruningMatrix getInvertedFromUnreducedPruningMatrix(SearchProblem sp) {
		ReducedPruningMatrix unreduced = new ReducedPruningMatrix(this);

		// invert the QStar pruneMat
		InvertedPruningMatrix ans = new InvertedPruningMatrix(sp, unreduced.getUpdatedPruningMatrix());

		// throw exception if there is a sequence mismatch
		if(!pruneMatIsValid(ans))
			throw new RuntimeException("ERROR: pruning did not reduce RCs to sequence space of allowedAAs");
		
		BigInteger numConfs = numConfs(ans);
		return numConfs.compareTo(BigInteger.ZERO) == 0 ? null : ans;
	}


	public InvertedPruningMatrix getInvertedFromReducedPruningMatrix(SearchProblem sp) {

		// invert the QStar pruneMat
		InvertedPruningMatrix ans = new InvertedPruningMatrix(sp,
				((ReducedPruningMatrix)sp.reducedMat).getUpdatedPruningMatrix());

		// throw exception if there is a sequence mismatch
		if(!pruneMatIsValid(ans))
			throw new RuntimeException("ERROR: pruning did not reduce RCs to sequence space of allowedAAs");

		BigInteger numConfs = numConfs(ans);
		return numConfs.compareTo(BigInteger.ZERO) == 0 ? null : ans;
	}


	public UnprunedPruningMatrix getUnprunedPruningMatrix( SearchProblem sp, double pruningInterval ) {

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
		for( int pos = 0; pos < sp.pruneMat.numPos(); ++pos ) {
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


	public SearchProblem getReducedSearchProblem( String name, 
			ArrayList<ArrayList<String>> allowedAAs, 
			ArrayList<String> flexRes, 
			ArrayList<Integer> posNums) {

		// Create a version of the search problem restricted to the specified sequence (list of amino acid names)
		// the constructor creates a new confspace object
		SearchProblem reducedSP = new SearchProblem(this, name, allowedAAs, flexRes, posNums);

		reducedSP.reducedMat = reducedSP.getReducedPruningMatrix(reducedSP); // for q*, q'
		
		if(reducedSP.reducedMat != null)
			reducedSP.inverseMat = reducedSP.getInvertedFromReducedPruningMatrix(reducedSP); // for p*

		return reducedSP;
	}


	private ArrayList<String> getAAsAtPos( PruningMatrix pruneMat, int pos ) {
		ArrayList<String> ans = new ArrayList<>();

		ArrayList<Integer> rcsAtPos = pruneMat.unprunedRCsAtPos(pos);
		for( int RCNum : rcsAtPos ) {
			int pos1 = posNums.get(pos);
			String AAType = confSpace.posFlex.get(pos1).RCs.get(RCNum).AAType;
			if(!ans.contains(AAType)) ans.add(AAType);
		}

		return ans;
	}


	private ArrayList<ArrayList<String>> getAAsAtPos( PruningMatrix pruneMat ) {
		ArrayList<ArrayList<String>> ans = new ArrayList<>();

		for( int pos = 0; pos < pruneMat.numPos(); ++pos )
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

		for( int pos = 0; pos < pruneMat.numPos(); ++pos ) {
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
}
