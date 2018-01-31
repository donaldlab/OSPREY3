/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupleEnumerator;
import edu.duke.cs.osprey.confspace.TupleMatrix;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.NewEPICMatrixCalculator;
import edu.duke.cs.osprey.ematrix.ReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.ematrix.epic.NewEPICMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.pruning.NewPruner;
import edu.duke.cs.osprey.pruning.NewPruningControl;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.BasicEPICTupleExpander;
import edu.duke.cs.osprey.tupexp.ConfETupleExpander;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import edu.duke.cs.osprey.tupexp.NewConfETupleExpander;
import edu.duke.cs.osprey.tupexp.TupExpChooser;
import edu.duke.cs.osprey.tupexp.TupleExpander;
import edu.duke.cs.osprey.voxq.VoxelGCalculator;

/**
 *
 * Information about RC tuples in a single-state design system (possibly part of a multistate design)
 * Stored in matrices.  Energy, pruning, etc. 
 * 
 * 
 * @author mhall44
 */
public class PrecomputedMatrices {
        
    
    String name;
    EnergyMatrix emat;//this may be handled separately?  For example it doesn't need to be updated
    PruningMatrix competitorPruneMat = null;
    PruningMatrix pruneMat = null;
    NewEPICMatrix epicMat = null;
    EnergyMatrix luteMat = null;
    
    SimpleConfSpace confSpace;
    EnergyCalculator ecalc;
    ConfEnergyCalculator confECalc;
    
    
    boolean EFullConfOnly = false;//for now!!
    
    
    EPICSettings epicSettings;
    LUTESettings luteSettings;
    PruningSettings pruningSettings;
    
    NewPruningControl pruningControl;
    double pruningInterval, Ew;

    public PrecomputedMatrices(double pruningInterval, double Ew, String name, EnergyMatrix emat, 
            SimpleConfSpace confSpace, EnergyCalculator ecalc, ConfEnergyCalculator confECalc,
            EPICSettings epicSettings, 
            LUTESettings luteSettings, PruningSettings pruningSettings) {
        this.name = name;
        this.emat = emat;
        this.confSpace = confSpace;
        this.confECalc = confECalc;
        this.ecalc = ecalc;
        this.epicSettings = epicSettings;
        this.luteSettings = luteSettings;
        this.pruningSettings = pruningSettings;
        this.pruningInterval = pruningInterval;
        this.Ew = Ew;
        precompute();
    }
    
    
    private void precompute(){
        if(EFullConfOnly){//Perform a tuple expansion that does not assume a pairwise
            //energy function, and thus must omit some of the usual pruning steps.
            fullConfOnlyTupExp();
            return;
        }
    
        //First calculate the pairwise energy matrix, if not already present
        if (emat == null) {
            loadEnergyMatrix();
        }
        
        
        //Set up pruning
        pruningControl = new NewPruningControl(
            this,
            0, // pruning interval, set by initPruning()
            pruningSettings.typedep,
            pruningSettings.boundsThresh,
            pruningSettings.algOption,
            pruningSettings.useFlags,
            pruningSettings.useTriples,
            false,
            false, // useEPIC, set by initPruning()
            false, // useTupExp, set by initPruning()
            pruningSettings.stericThresh
        );
        

        //Doing competitor pruning now
        //will limit us to a smaller, but effective, set of competitors in all future DEE
        if(competitorPruneMat == null){
            System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
            initPruning(0, false, false);
            //pruningControl.setOnlyGoldstein(true);//steric pruning essentially cuts conf space, so shouldn't
            //be competing with sterically pruned confs.  (Only an issue for questionable RCs of course)
            pruningControl.prune();
            competitorPruneMat = pruneMat;
            pruneMat = null;
            System.out.println("COMPETITOR PRUNING DONE");
        }
        
        
        //Next, do DEE, which will fill in the pruning matrix
        initPruning(pruningInterval, false, false);
        pruningControl.prune();//pass in DEE options, and run the specified types of DEE            
        
        
        //precomputing EPIC or tuple-expander matrices is much faster
        //if only done for unpruned RCs.  Less RCs to handle, and the fits are far simpler.  
        if(epicSettings.shouldWeUseEPIC()){
            loadEPICMatrix();
            
            //we can prune more using the EPIC matrix
            if(epicSettings.useEPICPruning){
                System.out.println("Beginning post-EPIC pruning.");
                initPruning(pruningInterval, true, false);
                pruningControl.prune();
                System.out.println("Finished post-EPIC pruning.");
            }
        }
        if(luteSettings.shouldWeUseLUTE()){//preferably do this one EPIC loaded (much faster if can fit to EPIC)
            loadTupExpEMatrix();
            
            //we can prune even more with tup-exp!
            //we can prune more using the EPIC matrix
            //no iMinDEE interval needed here
            System.out.println("Beginning post-tup-exp pruning.");
            initPruning(Ew, false, true);
            pruningControl.prune();
            System.out.println("Finished post-tup-exp pruning.");
        }
    }
    
 
    private void initPruning(double localPruningInterval, boolean useEPIC, boolean useTupExp) {
        
        // init the pruning matrix if needed
        if(pruneMat == null || pruneMat.getPruningInterval() < pruningInterval) {
            pruneMat = new PruningMatrix(confSpace.getNumPos(), confSpace.getNumResConfsByPos(), pruningInterval);
        }
        
        // configure the pruner
        pruningControl.setOnlyGoldstein(false);
        pruningControl.setPruningInterval(localPruningInterval);
        pruningControl.setUseEPIC(useEPIC);
        pruningControl.setUseTupExp(useTupExp);
    }
 
    
    private void fullConfOnlyTupExp(){
        //precompute the tuple expansion
        if(!luteSettings.shouldWeUseLUTE())
            throw new RuntimeException("ERROR: Need tuple expansion to handle full-conf-only E-function");
        if(epicSettings.shouldWeUseEPIC())//later consider using differencing scheme to do EPIC for these
            throw new RuntimeException("ERROR: EPIC for full-conf-only E-function not yet supported");
        
        
        //Let's compute a matrix from the pairwise terms (no P-B), to use in selecting triples
        loadEnergyMatrix();
        
        //initialize pruning matrix.  Nothing pruned yet because don't have pairwise energies
        pruneMat = new PruningMatrix(confSpace.getNumPos(), confSpace.getNumResConfsByPos(), pruningInterval);//not iMinDEE

        //We can only do steric pruning
        //May want to set a lower thresh than the default (30 perhaps)
        NewPruner pruner = new NewPruner(this, false, 0, 0, false, false);
        pruner.pruneSteric(pruningSettings.stericThresh);
                
        loadTupExpEMatrix();
    }
    
    
    public void loadEnergyMatrix(){
        loadMatrix(MatrixType.EMAT);
    }
    
    public void loadTupExpEMatrix(){
        loadMatrix(MatrixType.TUPEXPEMAT);
    }
    
    public void loadEPICMatrix(){
        loadMatrix(MatrixType.EPICMAT);
        
        //if(useVoxelG)
        //    gCalc = new VoxelGCalculator(this);
    }
    
    
    public enum MatrixType {
        EMAT, TUPEXPEMAT, EPICMAT;
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
            
        if(type == MatrixType.EMAT){
            return new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
        }
        else if(type == MatrixType.EPICMAT){
            NewEPICMatrixCalculator emCalc = new NewEPICMatrixCalculator(confSpace, confECalc, pruneMat, epicSettings);
            emCalc.calcPEM();
            return emCalc.getEPICMatrix();
            //EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(confSpace,shellResidues,pruneMat,epicSettings);
        }
        else {
            //need to calculate a tuple-expansion matrix
            
            //make a tuple expander
            TupleExpander expander;
            if(epicSettings.shouldWeUseEPIC()){
                expander = new BasicEPICTupleExpander(confSpace, pruningInterval, luteSettings,
                    epicMat, pruneMat);
            }
            else {
                expander = new NewConfETupleExpander(confSpace, pruningInterval, luteSettings,
                    confECalc, pruneMat);
            }
            
            TupleEnumerator tupEnum = new TupleEnumerator(pruneMat,emat,confSpace.getNumPos());
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
            epicMat = (NewEPICMatrix) matrixFromFile;
        else //tup-exp
            luteMat = (EnergyMatrix) matrixFromFile;
        
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

    public EnergyMatrix getEmat() {
        return emat;
    }

    public PruningMatrix getPruneMat() {
        return pruneMat;
    }

    public NewEPICMatrix getEpicMat() {
        return epicMat;
    }

    public EnergyMatrix getLuteMat() {
        return luteMat;
    }

    public PruningMatrix getCompetitorPruneMat() {
        return competitorPruneMat;
    }

    public SimpleConfSpace getConfSpace() {
        return confSpace;
    }
    
    public boolean shouldWeUseLUTE(){
        return luteSettings.shouldWeUseLUTE();
    }
    
    
    
    
    
}
