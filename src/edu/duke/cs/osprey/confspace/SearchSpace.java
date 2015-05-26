/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.util.ArrayList;
import java.util.Iterator;

/**
 *
 * @author mhall44
 */
public class SearchSpace {
    //This object keeps track of the positions and the possible assignments for them, as used in the search algorithms
    //generally these will be rDesidues (or super-residues) and their RCs; subclass SearchSpace to change this
    
    //We keep a ConfSpace together with "annotations" that help us find its GMEC, partition functions, etc.
    //annotations are based on RCs and indicate pairwise energies, pruning information, etc.
    //they also include iterators over RCs and pairs of interest
    
    public ConfSpace confSpace;
    
    public EnergyMatrix emat;//energy matrix.  Meanings:
    //-Defines full energy in the rigid or tuple-expander cases
    //-Defines ower bound in the classic continuous case
    //-emat + epicm = full energy if using EPIC for search 
    
    //EPICMatrix epicm = null;//EPIC matrix, to be used if appropriate
    //ADD SOON...
    
    public EnergyFunction fullConfE;//full energy for any conformation
    public ArrayList<Residue> shellResidues;//non-flexible residues to be accounted for in energy calculations
    
    public String name;//a human-readable name, which will also be used to name stored energy matrices, etc.
    
    public PruningMatrix pruneMat;
    
    boolean contSCFlex;
    /*
    //We'll also store iterators here
    //These iterators iterate over all RCs, and all pairs and higher tuples with nonnegligible energies
    Iterator<RC> RCIterator;
    Iterator<RC[]> pairIterator;
    Iterator<RC[]> higherIterator;
    
    //These iterators iterate only over nonpruned RCs and higher tuples
    Iterator<RC> unprunedRCIterator;//MAKE THIS A SPECIAL CLASS STORING LIST OF WHATS UNPRUNED
    //(for speed.  Used in A*.  May need one per position)
    
    RCHigherTupleIterator;
    *///probably not iterators just methods to get list of unpruned RCs, prunable tuples
    
    
    //some features to implement soon
    boolean useEPIC = false;
    boolean useTupExpForSearch = false;//use a tuple expansion to approximate the energy as we search
    
    boolean useEllipses = false; 
    
    public SearchSpace(
    		String name,
    		String PDBFile,
    		ArrayList<String> flexibleRes,
    		ArrayList<ArrayList<String>> allowedAAs,
    		boolean addWT,
            boolean contSCFlex,
            boolean ellipses){
    	
    	useEllipses = ellipses;
    	
        confSpace = new ConfSpace(PDBFile,flexibleRes,allowedAAs,addWT,contSCFlex, useEllipses);
        this.name = name;
        
        //initialize matrices
        emat = new EnergyMatrix(confSpace);
        pruneMat = new PruningMatrix(confSpace);//similar structure to emat but with booleans...
        
        //EPIC too!!
        
        this.contSCFlex = contSCFlex;
        
        
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
        return E;
    }
    
    public void outputMinimizedStruct(int[] conf, String PDBFileName){
        //Output minimized conformation to specified file
        //RCs are listed for all flexible positions in conf
        confSpace.minimizeEnergy(conf, fullConfE, PDBFileName);
    }
    
    
    public double approxMinimizedEnergy(int[] conf){
        //EPIC or other approximation for the minimized energy of the conformation
        //whose RCs are listed for all flexible positions in conf
        
        if(useEPIC){
            throw new RuntimeException("ERROR: EPIC not support yet...");
            /*double bound = emat.confE(conf);//emat contains the pairwise lower bounds
            double contPart = epicm.minContE(conf);//EPIC handles the continuous part (energy - pairwise lower bounds)
            
            return bound+contPart;
            */
        }
        else if( useTupExpForSearch || (!contSCFlex) ){//use matrix directly
            return emat.confE(conf);//matrix covers all energies
        }
        else
            throw new RuntimeException("ERROR: Asking searchSpace to approximate minimized energy but using a non-approximating method");
    }
    
    
    
    public double lowerBound(int[] conf){
        //Argument: RC assignments for all the flexible residues (RCs defined in resFlex)
        //return lower bound on energy for the conformational space defined by these assignments
        //based on precomputed energy matrix (including EPIC if indicated)
        
        double bound = emat.confE(conf);//the energy recorded by the matrix is 
        //the pairwise lower bounds
        
        return bound;
    }
    
    
    public void loadEnergyMatrix(){
        if(!loadEnergyMatrixFromFile()){
            emat = new EnergyMatrixCalculator(confSpace,shellResidues).calcPEM();
            storeEnergyMatrix();
            loadEnergyMatrixFromFile();
        }
    }
    
    
    boolean loadEnergyMatrixFromFile(){
        //try loading the pairwise energy matrix from a file
        //return true if successful, false if not, in which case we'll have to compute it
        //probably easiest to do this by serializing EnergyMatrix, possibly with custom serialization method...
        emat = (EnergyMatrix) ObjectIO.readObject(energyMatrixFileName(),true);
        return (emat!=null);//unsuccessful loading leaves null emat
    }
    
    
    void storeEnergyMatrix(){
        ObjectIO.writeObject(emat,energyMatrixFileName());
    }
    
    String energyMatrixFileName(){
        return name+".emat.dat";
    }
    
    //load EPICMatrix too...
    
    
    
}
