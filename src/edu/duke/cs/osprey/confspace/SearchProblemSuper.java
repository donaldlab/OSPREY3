/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
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
import java.util.Iterator;

/**
 *
 * @author hmn5
 */
public class SearchProblemSuper {
    //This is the super-RC analog of SearchProblem

    public ConfSpaceSuper confSpaceSuper;

    public EnergyMatrix emat;

    public EPICMatrix epicMat = null;
    public EPICSettings epicSettings = null;

    public EnergyMatrix tupExpEMat;

    public EnergyFunction fullConfE;//full energy for any conformation
    public ArrayList<Residue> shellResidues;

    public String name;

    public PruningMatrix pruneMat;

    boolean contSCFlex;

    public PruningMatrix competitorPruneMat;

    public boolean useEPIC = false;
    public boolean useTupExpForSearch = false;

    public SearchProblemSuper(SearchProblemSuper sp1) {
        //Shallow copy

        confSpaceSuper = sp1.confSpaceSuper;
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
    }

    public SearchProblemSuper(String name, String PDBFile, ArrayList<String> flexibleRes, ArrayList<ArrayList<String>> allowedAAs, boolean addWT,
            boolean contSCFlex, boolean useEPIC, EPICSettings epicSettings, boolean useTupExp, DEEPerSettings dset,
            ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, boolean useEllipses) {

        this.confSpaceSuper = new ConfSpaceSuper(PDBFile, flexibleRes, allowedAAs, addWT, contSCFlex, dset, moveableStrands, freeBBZones, useEllipses);

        this.name = name;

        this.contSCFlex = contSCFlex;
        this.useTupExpForSearch = useTupExp;
        this.useEPIC = useEPIC;
        this.epicSettings = epicSettings;

        //energy function setup
        EnergyFunctionGenerator eGen = EnvironmentVars.curEFcnGenerator;
        decideShellResidues(eGen.distCutoff);

        //fullConfE = eGen.fullConfEnergy(confSpace,shellResidues);
        //HMN: use confSpaceSuper as input instead of confSpace
        fullConfE = eGen.fullConfEnergy(confSpaceSuper, shellResidues);
    }

    private void decideShellResidues(double distCutoff) {
        //Decide what non-flexible residues need to be accounted for in energy calculations

        //use array list for these so they stay in order
        ArrayList<Residue> flexibleResidues = new ArrayList<>();
        for (PositionConfSpaceSuper pcs : confSpaceSuper.posFlex) {
            for (Residue res : pcs.resList) {
                flexibleResidues.add(res);
            }
        }

        //we'll decide what shell residues to include by using a simple distance cutoff with
        //the current conformation,
        //rather than doing a conformational search for the minimal distance (with respect to conformations)
        //of a shell residue to any flexible residues
        //the distance cutoff can be increased to accommodate this if desired.
        shellResidues = new ArrayList<>();

        for (Residue nonFlexRes : confSpaceSuper.m.residues) {
            //if residue is not flexible
            if (!flexibleResidues.contains(nonFlexRes)) {
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
        //Minimized eneryg of the conformatoin
        //whose super-RCs are listed for all flexible positions in conf
        double E = confSpaceSuper.minimizeEnergy(conf, fullConfE, null);
        return E;
    }

    public void outputMinimizedStruct(int[] conf, String PDBFileName) {
        //Output minimized conformation to specified file
        //super-RCs are listed for all flexible positions in conf
        confSpaceSuper.minimizeEnergy(conf, fullConfE, PDBFileName);
    }

    //HMN: Not currently supported
    public double approxMinimizedEnergy(int[] conf) {
        //EPIC or other approximation for the minimized energy of the conformation
        //whose RCs are listed for all flexible positions in conf
        throw new RuntimeException("TupExpansion and EPIC not currently supported with SearchProblemSuper and super-RCs");
        /*
         if (useTupExpForSearch) {//use tup-exp E-matrix direectly
         return tupExpEMat.confE(conf);
         } else if (useEPIC) {//EPIC w/o tup-exp
         return EPICMinimizedEnergy(conf);
         } else {
         throw new RuntimeException("ERROR: Asking searchSpace to approximate minimized energy but using a non-approximating method");
         }
         */
    }

    //HMN: Not currently supported
    public double EPICMinimizedEnergy(int[] conf) {
        throw new RuntimeException("TupExpansion and EPIC not currently supported with SearchProblemSuper and super-RCs");
        /*
         //approximate energy using EPIC
         double bound = emat.confE(conf);//emat contains the pairwise lower bounds
         double contPart = epicMat.minContE(conf);
         //EPIC handles the continuous part (energy - pairwise lower bounds)

         return bound + contPart;
         */
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
            EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpaceSuper, shellResidues);
            emCalc.calcPEM();
            return emCalc.getEMatrix();
        } else if (type == MatrixType.EPICMAT) {
            throw new RuntimeException("TupExpansion and EPIC not currently supported with SearchProblemSuper and super-RCs");
            /*
             EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpaceSuper, shellResidues,
             pruneMat, epicSettings);
             emCalc.calcPEM();
             return emCalc.getEPICMatrix();
             */
        } else {
            throw new RuntimeException("ONLY MatrixType.EMAT is currently supported in SearchProblemSuper and super-RCs");

            /*
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
             */
        }
    }

        boolean loadMatrixFromFile (MatrixType type, String matrixFileName) {
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
    }
