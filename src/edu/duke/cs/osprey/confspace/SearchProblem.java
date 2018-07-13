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

package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.util.ArrayList;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.ReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.ConfETupleExpander;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import edu.duke.cs.osprey.tupexp.TupExpChooser;
import edu.duke.cs.osprey.voxq.VoxelGCalculator;

/**
 *
 * @author mhall44
 */
public class SearchProblem implements Serializable {
    //This object keeps track of the positions and the possible assignments for them, as used in the search algorithms
    //generally these will be residues (or super-residues) and their RCs; subclass SearchProblem to change this
    
    //We keep a ConfSpace together with "annotations" that help us find its GMEC, partition functions, etc.
    //annotations are based on RCs and indicate pairwise energies, pruning information, etc.
    //they also include iterators over RCs and pairs of interest
    
	private static final long serialVersionUID = 2590525329048496524L;

    
    public ConfSpace confSpace;
    
    public EnergyMatrix emat;//energy matrix.  Meanings:
    //-Defines full energy in the rigid case
    //-Defines lower bound in the continuous case
    //-emat + epicm = full energy if using EPIC for search 
    
    public EPICMatrix epicMat = null;//EPIC matrix, to be used if appropriate
    public EPICSettings epicSettings = null;
    public LUTESettings luteSettings = null;
    public DEEPerSettings deeperSettings = null;
    
    public EnergyMatrix tupExpEMat;//Defines full energy in the continuous, tuple-expander case
    
    public EnergyFunction fullConfE;//full energy for any conformation
    public ArrayList<Residue> shellResidues;//non-flexible residues to be accounted for in energy calculations
    
    public String name;//a human-readable name, which will also be used to name stored energy matrices, etc.
    
    public PruningMatrix pruneMat;
    
    public boolean contSCFlex;
    public boolean useVoxelG = false;//use the free energy of each voxel instead of minimized energy
    VoxelGCalculator gCalc = null;
    
    public PruningMatrix competitorPruneMat;//a pruning matrix performed at pruning interval 0,
    //to decide which RC tuples are valid competitors for pruning
    
    public ArrayList<String> flexRes; 
    public ArrayList<ArrayList<String>> allowedAAs;
    
    public boolean useEPIC = false;
    public boolean useTupExpForSearch = false;//use a tuple expansion to approximate the energy as we search
    
    
    public boolean useERef = false;
    public boolean addResEntropy = false;
    
    public int numEmatThreads = 1;

    public PolytopeMatrix plugMat;

    
    public SearchProblem(SearchProblem other){//shallow copy
    	confSpace = other.confSpace;
    	emat = other.emat;
        epicMat = other.epicMat;
        epicSettings = other.epicSettings;
        luteSettings = other.luteSettings;
        deeperSettings = other.deeperSettings;
        tupExpEMat = other.tupExpEMat;
        
    	fullConfE = other.fullConfE;
    	shellResidues = other.shellResidues;
    	name = other.name + System.currentTimeMillis();//probably will want to change this to something more meaningful
        
    	pruneMat = other.pruneMat;
    	competitorPruneMat = other.competitorPruneMat;
        
    	contSCFlex = other.contSCFlex;
    	flexRes = other.flexRes;
        allowedAAs = other.allowedAAs;
    	
    	useVoxelG = other.useVoxelG;
        gCalc = other.gCalc;
    	
    	useEPIC = other.useEPIC;
    	useTupExpForSearch = other.useTupExpForSearch;
        
        useERef = other.useERef;
        addResEntropy = other.addResEntropy;
        numEmatThreads = other.numEmatThreads;
    }
    
    
    
    public SearchProblem(String name, String PDBFile, ArrayList<String> flexibleRes, ArrayList<ArrayList<String>> allowedAAs, boolean addWT,
            boolean contSCFlex, boolean useEPIC, EPICSettings epicSettings, boolean useTupExp, LUTESettings luteSettings, DEEPerSettings dset, 
            ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, boolean useEllipses, boolean useERef,
            boolean addResEntropy, boolean addWTRots, ResidueTermini termini, boolean useVoxelG, ArrayList<String> wtRotOnlyRes){
        
        confSpace = new ConfSpace(PDBFile, flexibleRes, allowedAAs, addWT, wtRotOnlyRes,
                contSCFlex, dset, moveableStrands, freeBBZones, useEllipses, addWTRots, termini);
        this.name = name;
        this.flexRes = flexibleRes;
        this.allowedAAs = allowedAAs;
        
        this.contSCFlex = contSCFlex;
        this.useTupExpForSearch = useTupExp;
        this.useEPIC = useEPIC;
        this.epicSettings = epicSettings;
        this.luteSettings = luteSettings;
        this.deeperSettings = dset;
        
        this.useERef = useERef;
        this.addResEntropy = addResEntropy;
        this.useVoxelG = useVoxelG;
        
        //energy function setup
        EnergyFunctionGenerator eGen = EnvironmentVars.curEFcnGenerator;
        decideShellResidues(eGen.ffParams.shellDistCutoff);
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
        if(useVoxelG)//use free instead of minimized energy
            return gCalc.calcG(conf);
        
        double E = confSpace.minimizeEnergy(conf, fullConfE, null);
        
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
        
        if( useTupExpForSearch ){//use tup-exp E-matrix directly
            return tupExpEMat.confE(conf);
        }
        else if( useEPIC ){//EPIC w/o tup-exp
            return EPICMinimizedEnergy(conf);
        }
        
        else
            throw new RuntimeException("ERROR: Asking searchSpace to approximate minimized energy but using a non-approximating method");
    }
    
    
    public double voxelFreeEnergy(int[] conf){
        if(gCalc==null)
            throw new RuntimeException("ERROR: Free energy calculator is null (probably no EPIC matrix loaded)");
            
        return gCalc.calcG(conf);
    }
    
    public double EPICMinimizedEnergy(int[] conf){
        //approximate energy using EPIC
        if(useVoxelG)
            return voxelFreeEnergy(conf);
        
        return epicMat.minimizeEnergy(new RCTuple(conf), true);
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
        
        if(useVoxelG)
            gCalc = new VoxelGCalculator(this);
    }

    public void loadPLUGMatrix(){//just loads the matrix--pruning can be handled by multi-term pruner
        loadMatrix(MatrixType.PLUGMAT);
    }
    
    
    public enum MatrixType {
        EMAT, TUPEXPEMAT, EPICMAT, PLUGMAT;
    }
    
    
    //load the specified matrix; if the right file isn't available then compute and store it
    public void loadMatrix(MatrixType type){
        
        String matrixFileName = name + "." + type.name() + ".dat";
        //matrix file names are determined by the name of the search problem
        
        if(!loadMatrixFromFile( type, matrixFileName )){
            TupleMatrix<?> matrix = calcMatrix(type);
            ObjectIO.writeObject( matrix, matrixFileName );
            loadMatrixFromFile( type, matrixFileName );
        }
    }
    
    
    //compute the matrix of the specified type
    public TupleMatrix<?> calcMatrix(MatrixType type){
    
        // TODO: the search problem shouldn't concern itself with energy matrices and how to compute them
        
        if(type == MatrixType.EMAT){
            //MH: All the current conformational perturbations as of 9/12/16 should support copying
            //to new molecules, but I'll leave this option in case new DOFs or something cause an issue
        	            
            // if we're using MPI or DEEPer, use the old energy matrix calculator
            if (EnvironmentVars.useMPI || (deeperSettings != null && deeperSettings.doPerturbations())) {
            	System.out.println("\n\nWARNING: concurrent minimizations disabled\n"); 
            	
                // see if the user tried to use threads too for some reason and try to be helpful
                if (EnvironmentVars.useMPI && numEmatThreads > 1) {
                    System.out.println("\n\nWARNING: multiple threads and MPI both configured for emat calculation."
                        + " Ignoring thread settings and using only MPI.\n");
                }
                
                EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpace, shellResidues, useERef, addResEntropy);
                emCalc.calcPEM();
                return emCalc.getEMatrix();
            
            } else {
            
                // otherwise, use the new multi-threaded calculator (which doesn't support MPI)
                
                // calculate the emat! Yeah!
                SimpleEnergyMatrixCalculator ecalc = new SimpleEnergyMatrixCalculator.Cpu(numEmatThreads, EnvironmentVars.curEFcnGenerator.ffParams, confSpace, shellResidues);
                EnergyMatrix emat = ecalc.calcEnergyMatrix();
                
                // need to subtract reference energies?
                if (useERef) {
                    System.out.println("Computing reference energies...");
                    emat.seteRefMat(new ReferenceEnergies(confSpace));
                }
                
                // need to add entropies?
                if (addResEntropy) {
                    System.out.println("Computing residue entropies...");
                	for (int pos=0; pos<emat.getNumPos(); pos++) {
                		for (int rc=0; rc<emat.getNumConfAtPos(pos); rc++) {
                			double energy = emat.getOneBody(pos, rc);
                			energy += confSpace.getRCResEntropy(pos, rc);
                			emat.setOneBody(pos, rc, energy);
                		}
                	}
                }
                
                // cleanup
                ecalc.cleanup();
                
                return emat;
            }
        }
        else if(type == MatrixType.EPICMAT){
            EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpace,shellResidues,
                    pruneMat,epicSettings);
            emCalc.calcPEM();
            return emCalc.getEPICMatrix();
        }
        else if(type == MatrixType.PLUGMAT){//Let's not prune for now
            return new PolytopeMatrix(this, false);
        }
        else {
            //need to calculate a tuple-expansion matrix
            
            ConfETupleExpander expander = new ConfETupleExpander(this);//make a tuple expander
            
            
            TupleEnumerator tupEnum = new TupleEnumerator(pruneMat,emat,confSpace.numPos);
            TupExpChooser chooser = new TupExpChooser(expander, tupEnum);//make a chooser to choose what tuples will be in the expansion
            
            double curResid = chooser.calcPairwiseExpansion();//start simple...
            
            if(curResid > luteSettings.goalResid){//go to triples if needed
                System.out.println("EXPANDING PAIRWISE EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (2 PARTNERS)...");
                curResid = chooser.calcExpansionResTriples(2);
            }
            if(curResid > luteSettings.goalResid){//go to 5 partners if still need better resid...
                System.out.println("EXPANDING EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (5 PARTNERS)...");
                curResid = chooser.calcExpansionResTriples(5);
            }
            if(curResid > luteSettings.goalResid){
                System.out.println("WARNING: Desired LUTE residual threshold "+
                        luteSettings.goalResid+" not reached; best="+curResid);
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
        else if(type == MatrixType.PLUGMAT)
            plugMat = (PolytopeMatrix) matrixFromFile;
        else //tup-exp
            tupExpEMat = (EnergyMatrix) matrixFromFile;
        
        if(matrixFromFile==null)//unsuccessful loading leaves null emat
            return false;
        
        
        //check pruning interval.  Current interval is in pruneMat if we have pruned already;
        //if not then we need a matrix with infinite pruning interval (valid for all RCs).
        double matrixPruningInterval = ((TupleMatrix<?>)matrixFromFile).getPruningInterval();
        
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
    
    
    
    public boolean searchNeedsHigherOrderTerms(){
        //Will a conformation search using this SearchProblem need to use
        //higher-than-pairwise terms?
        if(useTupExpForSearch)
            return tupExpEMat.hasHigherOrderTerms();
        else
            return emat.hasHigherOrderTerms();
    }
   
}
