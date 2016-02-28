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
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.ConfETupleExpander;
import edu.duke.cs.osprey.tupexp.TupExpChooser;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
@SuppressWarnings("serial")
public class SearchProblem implements Serializable {
    //This object keeps track of the positions and the possible assignments for them, as used in the search algorithms
    //generally these will be rDesidues (or super-residues) and their RCs; subclass SearchProblem to change this
    
    //We keep a ConfSpace together with "annotations" that help us find its GMEC, partition functions, etc.
    //annotations are based on RCs and indicate pairwise energies, pruning information, etc.
    //they also include iterators over RCs and pairs of interest
    
    public ConfSpace confSpace;
    
    public EnergyMatrix emat;//energy matrix.  Meanings:
    //-Defines full energy in the rigid case
    //-Defines lower bound in the continuous case
    //-emat + epicm = full energy if using EPIC for search 
    
    public EPICMatrix epicMat = null;//EPIC matrix, to be used if appropriate
    public EPICSettings epicSettings = null;
    
    public EnergyMatrix tupExpEMat;//Defines full energy in the continuous, tuple-expander case
    
    public EnergyFunction fullConfE;//full energy for any conformation
    public ArrayList<Residue> shellResidues;//non-flexible residues to be accounted for in energy calculations
    
    public String name;//a human-readable name, which will also be used to name stored energy matrices, etc.
    
    public PruningMatrix pruneMat;
    
    boolean contSCFlex;
    
    public ArrayList<String> flexibleRes;
	public ArrayList<ArrayList<String>> allowedAAs;
	public String PDBFile;
    
    public PruningMatrix competitorPruneMat;//a pruning matrix performed at pruning interval 0,
    //to decide which RC tuples are valid competitors for pruning
    
        
    
    public boolean useEPIC = false;
    public boolean useTupExpForSearch = false;//use a tuple expansion to approximate the energy as we search
    
    
    boolean useERef = false;
    boolean addResEntropy = false;
    
    
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
    }
    
    
    
    public SearchProblem(String name, String PDBFile, ArrayList<String> flexibleRes, ArrayList<ArrayList<String>> allowedAAs, boolean addWT,
            boolean contSCFlex, boolean useEPIC, EPICSettings epicSettings, boolean useTupExp, DEEPerSettings dset, 
            ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, boolean useEllipses, boolean useERef,
            boolean addResEntropy){
        
        confSpace = new ConfSpace(PDBFile, flexibleRes, allowedAAs, addWT, contSCFlex, dset, moveableStrands, freeBBZones, useEllipses);
        this.name = name;
        
        this.flexibleRes = flexibleRes;
		this.allowedAAs = allowedAAs;
		this.PDBFile = PDBFile;
        
        this.contSCFlex = contSCFlex;
        this.useTupExpForSearch = useTupExp;
        this.useEPIC = useEPIC;
        this.epicSettings = epicSettings;
        
        this.useERef = useERef;
        this.addResEntropy = addResEntropy;
        
        //energy function setup
        EnergyFunctionGenerator eGen = EnvironmentVars.curEFcnGenerator;
        decideShellResidues(eGen.distCutoff);
        fullConfE = eGen.fullConfEnergy(confSpace,shellResidues);
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
        
        if(useERef)
            E -= emat.geteRefMat().confERef(conf);
        
        if(addResEntropy)
            E += confSpace.getConfResEntropy(conf);            
        
        return E;
    }
    
    
    public double rigidEnergy(int[] conf) {
    	//Minimized energy of the conformation
        //whose RCs are listed for all flexible positions in conf
        double E = confSpace.rigidEnergy(conf, fullConfE);
        
        if(useERef)
            E -= emat.geteRefMat().confERef(conf);
        
        if(addResEntropy)
            E += confSpace.getConfResEntropy(conf);            
        
        return E;
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
    
    
    
    
    
    //LOADING AND PRECOMPUTATION OF ENERGY MATRIX-TYPE OBJECTS (regular energy matrix, tup-exp and EPIC matrices)
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
    
	public ArrayList<String> getFlexibleResiduePositions(ArrayList<String> seq, ArrayList<Integer> ordinalPos){
		// converts ordinal position to absolute position in the molecule
		ArrayList<String> absolutePos = new ArrayList<>();

		for( int i = 0; i < ordinalPos.size(); ++i ) {
			int pos = ordinalPos.get(i);
			absolutePos.add(this.flexibleRes.get(pos));

			String aa = seq.get(i);
			if( !this.allowedAAs.get(pos).contains(aa) )
				throw new RuntimeException("ERROR: the specified amino acid " + aa 
						+ " at oridinal position " + i + " is not allowed");
		}

		return absolutePos;
	}
    
}
