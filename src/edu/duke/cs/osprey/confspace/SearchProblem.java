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
import java.util.Collections;
import java.util.HashMap;

/**
 *
 * @author mhall44
 */
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

    public boolean contSCFlex;

    public PruningMatrix competitorPruneMat;//a pruning matrix performed at pruning interval 0,
    //to decide which RC tuples are valid competitors for pruning

    public boolean useEPIC = false;
    public boolean useTupExpForSearch = false;//use a tuple expansion to approximate the energy as we search

    boolean useERef = false;

    public SearchProblem(SearchProblem sp1) {//shallow copy
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

    }

    public SearchProblem() {
    }

    public SearchProblem(String name, String PDBFile, ArrayList<String> flexibleRes, ArrayList<ArrayList<String>> allowedAAs, boolean addWT,
            boolean addWTRots, boolean contSCFlex, boolean useEPIC, EPICSettings epicSettings, boolean useTupExp, DEEPerSettings dset,
            ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, boolean useEllipses, boolean useERef) {

        confSpace = new ConfSpace(PDBFile, flexibleRes, allowedAAs, addWT, addWTRots, contSCFlex, dset, moveableStrands, freeBBZones, useEllipses);
        this.name = name;

        this.contSCFlex = contSCFlex;
        this.useTupExpForSearch = useTupExp;
        this.useEPIC = useEPIC;
        this.epicSettings = epicSettings;

        this.useERef = useERef;
        //energy function setup
        EnergyFunctionGenerator eGen = EnvironmentVars.curEFcnGenerator;
        decideShellResidues(eGen.distCutoff);
        fullConfE = eGen.fullConfEnergy(confSpace, shellResidues);
    }

    private void decideShellResidues(double distCutoff) {
        //Decide what non-flexible residues need to be accounted for in energy calculations

        ArrayList<Residue> flexibleResidues = new ArrayList<>();//array list for these so they stay in order
        for (PositionConfSpace pcs : confSpace.posFlex) {
            flexibleResidues.add(pcs.res);
        }

        //we'll decide what shell residues to include by using a simple distance cutoff with
        //the current conformation,
        //rather than doing a conformational search for the minimal distance (with respect to conformations)
        //of a shell residue to any flexible residues
        //the distance cutoff can be increased to accommodate this if desired.
        shellResidues = new ArrayList<>();

        for (Residue nonFlexRes : confSpace.m.residues) {
            if (!flexibleResidues.contains(nonFlexRes)) {//residue is not flexible

                for (Residue flexRes : flexibleResidues) {
                    double dist = flexRes.distanceTo(nonFlexRes);
                    if (dist <= distCutoff) {
                        shellResidues.add(nonFlexRes);//close enough to flexRes that we should add it
                        break;
                    }
                }
            }
        }
    }

    public double minimizedEnergy(int[] conf) {
        //Minimized energy of the conformation
        //whose RCs are listed for all flexible positions in conf
        double E = confSpace.minimizeEnergy(conf, fullConfE, null);

        if (useERef) {
            E -= emat.geteRefMat().confERef(conf);
        }

        return E;
    }

    public void outputMinimizedStruct(int[] conf, String PDBFileName) {
        //Output minimized conformation to specified file
        //RCs are listed for all flexible positions in conf
        //Note: eRef not included (just minimizing w/i voxel)
        confSpace.minimizeEnergy(conf, fullConfE, PDBFileName);
    }

    public double approxMinimizedEnergy(int[] conf) {
        //EPIC or other approximation for the minimized energy of the conformation
        //whose RCs are listed for all flexible positions in conf

        if (useTupExpForSearch) {//use tup-exp E-matrix direectly
            return tupExpEMat.confE(conf);
        } else if (useEPIC) {//EPIC w/o tup-exp
            return EPICMinimizedEnergy(conf);
        } else {
            throw new RuntimeException("ERROR: Asking searchSpace to approximate minimized energy but using a non-approximating method");
        }
    }

    public double EPICMinimizedEnergy(int[] conf) {
        //approximate energy using EPIC
        double bound = emat.confE(conf);//emat contains the pairwise lower bounds
        double contPart = epicMat.minContE(conf);
        //EPIC handles the continuous part (energy - pairwise lower bounds)

        return bound + contPart;
    }

    public double lowerBound(int[] conf) {
        //Argument: RC assignments for all the flexible residues (RCs defined in resFlex)
        //return lower bound on energy for the conformational space defined by these assignments
        //based on precomputed energy matrix (including EPIC if indicated)

        double bound = emat.confE(conf);//the energy recorded by the matrix is 
        //the pairwise lower bounds

        return bound;
    }

    //LOADING AND PRECOMPUTATION OF ENERGY MATRIX-TYPE OBJECTS (regular energy matrix, tup-exp and EPIC matrices)
    public void loadEnergyMatrix() {
        loadMatrix(MatrixType.EMAT);
    }

    public void loadTupExpEMatrix() {
        loadMatrix(MatrixType.TUPEXPEMAT);
    }

    public void loadEPICMatrix() {
        loadMatrix(MatrixType.EPICMAT);
    }

    enum MatrixType {

        EMAT, TUPEXPEMAT, EPICMAT;
    }

    //load the specified matrix; if the right file isn't available then compute and store it
    private void loadMatrix(MatrixType type) {

        String matrixFileName = name + "." + type.name() + ".dat";
        //matrix file names are determined by the name of the search problem

        if (!loadMatrixFromFile(type, matrixFileName)) {
            TupleMatrix matrix = calcMatrix(type);
            ObjectIO.writeObject(matrix, matrixFileName);
            loadMatrixFromFile(type, matrixFileName);
        }
    }

    //compute the matrix of the specified type
    private TupleMatrix calcMatrix(MatrixType type) {

        if (type == MatrixType.EMAT) {
            EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpace, shellResidues, useERef);
            emCalc.calcPEM();
            return emCalc.getEMatrix();
        } else if (type == MatrixType.EPICMAT) {
            EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpace, shellResidues,
                    pruneMat, epicSettings);
            emCalc.calcPEM();
            return emCalc.getEPICMatrix();
        } else {
            //need to calculate a tuple-expansion matrix

            double errorThresh = 0.01;

            ConfETupleExpander expander = new ConfETupleExpander(this);//make a tuple expander
            TupleEnumerator tupEnum = new TupleEnumerator(pruneMat, emat, confSpace.numPos);
            TupExpChooser chooser = new TupExpChooser(expander, tupEnum);//make a chooser to choose what tuples will be in the expansion

            double curResid = chooser.calcPairwiseExpansion();//start simple...

            if (curResid > errorThresh) {//go to triples if needed
                System.out.println("EXPANDING PAIRWISE EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (2 PARTNERS)...");
                curResid = chooser.calcExpansionResTriples(2);
            }
            if (curResid > errorThresh) {//go to 5 partners if still need better resid...
                System.out.println("EXPANDING EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (5 PARTNERS)...");
                curResid = chooser.calcExpansionResTriples(5);
            }
            if (curResid > errorThresh) {
                System.out.println("WARNING: Desired LUTE residual threshold "
                        + errorThresh + " not reached; best=" + curResid);
            }

            return expander.getEnergyMatrix();//get the final energy matrix from the chosen expansion
        }
    }

    boolean loadMatrixFromFile(MatrixType type, String matrixFileName) {
        //try loading the specified matrix from a file
        //return true if successful, false if not, in which case we'll have to compute it
        //also if the matrix's pruning interval is too low, it may be missing some RCs
        //that are unpruned at our current pruningInterval, so we have to recompute
        Object matrixFromFile = ObjectIO.readObject(matrixFileName, true);

        if (type == MatrixType.EMAT) {
            emat = (EnergyMatrix) matrixFromFile;
        } else if (type == MatrixType.EPICMAT) {
            epicMat = (EPICMatrix) matrixFromFile;
        } else //tup-exp
        {
            tupExpEMat = (EnergyMatrix) matrixFromFile;
        }

        if (matrixFromFile == null)//unsuccessful loading leaves null emat
        {
            return false;
        }

        //check pruning interval.  Current interval is in pruneMat if we have pruned already;
        //if not then we need a matrix with infinite pruning interval (valid for all RCs).
        double matrixPruningInterval = ((TupleMatrix) matrixFromFile).getPruningInterval();

        if (matrixPruningInterval == Double.POSITIVE_INFINITY)//definitely valid
        {
            return true;
        } else {
            //excludes some RCs...check against pruneMat pruning interval
            if (pruneMat == null) {
                throw new RuntimeException("ERROR: Trying to load pruning-dependent tuple matrix"
                        + "(EPIC or tup-exp) but haven't pruned yet");
            }

            return (matrixPruningInterval >= pruneMat.getPruningInterval());
        }
    }

    /**
     * HMN: Returns a search problem over a partial search space defined by a
     * subset of the total number of positions
     *
     * @param subsetOfPositions the subset of positions
     * @param updatedPruneMat updated pruning matrix (either empty or has one prune mat)
     * @return
     */
    public SearchProblem getPartialSearchProblem(ArrayList<Integer> subsetOfPositions, PruningMatrix... updatePruneMat) {
        //Make sure the positions are sorted
        Collections.sort(subsetOfPositions);

        //Create a shallow copy 
        SearchProblem partSearchSpace = new SearchProblem(this);
        //Update the pruning matrix
        if (updatePruneMat.length != 0){
            partSearchSpace.pruneMat = updatePruneMat[0];
        }
        
        ConfSpace partCSpace = partSearchSpace.confSpace.getPartialConfSpace(subsetOfPositions);
        EnergyMatrix partEmat = new EnergyMatrix(partSearchSpace.emat.getSubsetMatrix(subsetOfPositions));
        PruningMatrix partPruneMat = new PruningMatrix(partSearchSpace.pruneMat.getSubsetMatrix(subsetOfPositions));
        PruningMatrix partCompetitorPruneMat = new PruningMatrix(partSearchSpace.competitorPruneMat.getSubsetMatrix(subsetOfPositions));

        EPICMatrix partEPICMat = null;
        EnergyMatrix partTupExpEMat = null;
        if (useEPIC) {
            partEPICMat = new EPICMatrix(partSearchSpace.epicMat.getSubsetMatrix(subsetOfPositions), partCSpace);
        }
        if (useTupExpForSearch) {
            partTupExpEMat = new EnergyMatrix(partSearchSpace.tupExpEMat.getSubsetMatrix(subsetOfPositions));
        }
        //Currently, the partial confSpace contains a new molecule since it was 
        //deep copied. Thus, we need to get the shellResidues that correspond to
        //this new molecules
        //We will add all shellResidues from the full Search Problem (this could
        //change if we want to incorporate residue specific energy-cut offs)
        ArrayList<Residue> partShellResidues = new ArrayList<Residue>();
        //Create a hashmap for efficient adding
        //Maps resNum to the Residue
        HashMap<Integer, Residue> resNumToResidue = new HashMap<>();
        //For every residue in the molecules we add the information to the hashmap
        partCSpace.m.residues.stream().map(res -> resNumToResidue.put(res.resNum, res));
        //For each original shell residue, add the corresponding residue from
        //our partCSpace molecule to the list
        this.shellResidues.stream().map(res -> partShellResidues.add(resNumToResidue.get(res.resNum)));

        EnergyFunction partFullConfE = EnvironmentVars.curEFcnGenerator.fullConfEnergy(partCSpace, partShellResidues);

        //Now lets update partSearchSpace with the new information
        partSearchSpace.confSpace = partCSpace;
        partSearchSpace.emat = partEmat;
        partSearchSpace.epicMat = partEPICMat;
        partSearchSpace.tupExpEMat = partTupExpEMat;
//        partSearchSpace.name = aName;
        partSearchSpace.name = this.name + System.currentTimeMillis();
        partSearchSpace.pruneMat = partPruneMat;
        partSearchSpace.competitorPruneMat = partCompetitorPruneMat;

        return partSearchSpace;
    }

    /**
     * HMN: This is used for partial search spaces. Given an interaction graph,
     * we update the underlying "energy" matrices to only have non-zero entries
     * for the pairwise terms defined by the interaction graph
     *
     * @param interactionGraph the interaction graph that maps position_I and
     * position_J to true if we should keep the pairwise energies between
     * position_I and position_J
     */
    public void updateMatrixCrossTerm(boolean[][] interactionGraph) {
        this.emat.updateMatrixCrossTerms(interactionGraph);
        if (useEPIC) {
            this.epicMat.updateMatrixCrossTerms(interactionGraph);
        }
        if (useTupExpForSearch) {
            this.tupExpEMat.updateMatrixCrossTerms(interactionGraph);
        }
    }

    
    /**
     * HMN: This method is subtracts the internal energies from another matrix, 
     * presumably the unbound matrix for partial search spaces
     */
    public void substractUnboundInternalEnergies(SearchProblem unBoundSearchProblem, ArrayList<Integer> posNumsToSubtractFrom, HashMap<Integer, Integer> boundPosNumToUnboundPosNum) {
        this.emat.subtractUnboundInternalEnergies(unBoundSearchProblem.emat, posNumsToSubtractFrom, boundPosNumToUnboundPosNum);
        if (contSCFlex){
            throw new RuntimeException("Continuous Flexibility Not Fully Implemented in KaDEE");
        }
        if (useEPIC) {
            throw new RuntimeException("Subtracting Unbound Energies Currently not supported with EPIC");
        }
        if (useTupExpForSearch) {
            this.tupExpEMat.subtractUnboundInternalEnergies(unBoundSearchProblem.tupExpEMat, posNumsToSubtractFrom, boundPosNumToUnboundPosNum);
        }
    }
    
}
