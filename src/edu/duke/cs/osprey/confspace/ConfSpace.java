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
import java.math.BigInteger;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.MoveableStrand;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.dof.StrandRotation;
import edu.duke.cs.osprey.dof.StrandTranslation;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.StringParsing;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author mhall44
 */
public class ConfSpace implements Serializable {
    //This class represents the conformational search space for a design
    //used for GMEC-based design, K*, or anything else we want to do with a conformational space
    //This class just defines the conformational space itself (the molecule + all kinds of flexibility
    //and possible mutations, and how these are represented as RCs, etc.)
    //This class can be put in a SearchProblem to add annotations like what RCs are pruned,
    //what their pairwise energies are, etc.  
    
	private static final long serialVersionUID = 6414329117813457771L;    
 
	public Molecule m;
    //The molecule will be composed of residues. 
    //It will have one set of coordinates, which are stored in the residues to make mutation
    //and pairwise energy computation easy (no need for complicated partial arrays, subtracting off
    //template energies, etc., and mutation will be very quick and will only require
    //AMBER reinitialization for the affected residue)
        
    //If we need to keep a copy of, say, the original coordinates, we can have a second molecule origMolec
    //for that
    //Once loaded, the molecule can only be changed by functions overriding DegreeOfFreedom.applyValue
    //this will give us a uniform framework for applying conformational and sequence changes
    //and each DOF will store its value, so we can easily and unambiguously look up the state of the
    //molecule at any given time
    
    
    public ArrayList<DegreeOfFreedom> confDOFs = new ArrayList<>();
    //list of conformational degrees of freedom
    
    public ArrayList<ResidueTypeDOF> mutDOFs = new ArrayList<>();
    //list of sequence degrees of freedom (residue types for the mutable residues)
    
    
    public ArrayList<PositionConfSpace> posFlex = new ArrayList<>();
    //defines the flexible positions and what RCs they have
    //generally each position is a residue, but it could be more than one ("super-residue" with "super-RCs")

    public ArrayList<String> flexibleRes;
    public int numPos;//number of flexible positions
    
    public boolean useEllipses = false;
    
    /** initialize a new conformational space, defining all its flexibility
     *   we use one residue per position here
     * 
     * @param PDBFile the structure to read from
     * @param flexibleRes list of residue numbers to be made flexible (as in PDB file)
     * @param allowedAAs list of allowed residue types at each flexible position
     * @param addWT whether to add wild-type to the allowed AA types 
     * @param contSCFlex means allow continuous sidechain flexibility
     * @param dset DEEPer Settings
     * @param moveableStrands ... ? 
     * @param freeBBZones ...? 
     * @param ellipses model ellipses
     * @param addWTRots add the wild-type 'rotamers'
     */
    public ConfSpace(String PDBFile, ArrayList<String> flexibleRes, ArrayList<ArrayList<String>> allowedAAs, 
            boolean addWT, ArrayList<String> wtRotOnlyRes, boolean contSCFlex, DEEPerSettings dset, ArrayList<String[]> moveableStrands, 
            ArrayList<String[]> freeBBZones, boolean ellipses, boolean addWTRots, ResidueTermini termini){
    
    	useEllipses = ellipses;  	
    	this.flexibleRes = flexibleRes;
        numPos = flexibleRes.size();
        
        //read the structure and assign templates, deleting unassignable res...
        m = new Strand.Builder(PDBIO.readFile(PDBFile)).setResidues(termini).build().mol;
        
        // before making any structure changes, capture the wt rots if needed
        List<ResidueTemplate> wtRots = new ArrayList<>(Collections.nCopies(numPos, null));
        if (addWTRots) {
        	for (int i=0; i<numPos; i++) {
        		// TODO: support alternate conformations?
        		Residue res = m.getResByPDBResNumber(flexibleRes.get(i));
        		wtRots.set(i, ResidueTemplate.makeFromResidueConfs(res));
        	}
        }
        
        //Make all the degrees of freedom
        //start with proline puckers (added to res)
        makeRingPuckers(allowedAAs, flexibleRes);
        
        //now the single-res degrees of freedom (mutations, sidechain dihedrals)
        ArrayList<ArrayList<DegreeOfFreedom>> singleResDOFs = new ArrayList<>();

        for(int pos=0; pos<numPos; pos++){
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            
            //at this point, m has all wild-type residues, so just see what res is now
            String wtName = res.template.name;
            if( ! StringParsing.containsIgnoreCase(allowedAAs.get(pos), wtName) ){//wtName not currently in allowedAAs
                if(addWT || allowedAAs.get(pos).isEmpty())//It should be (addWT or blank AA type)
                    allowedAAs.get(pos).add(wtName);
                else//it should not be...make sure wt rots aren't included
                    wtRots.set(pos, null);
            }
            
            ArrayList<DegreeOfFreedom> resDOFs = mutablePosDOFs(res,allowedAAs.get(pos));//add mutation and dihedral confDOFs
                        
            ResidueTypeDOF resMutDOF = (ResidueTypeDOF)resDOFs.remove(0);//first mutable pos DOF is the mutation-type DOF
 
            mutDOFs.add(resMutDOF);

            singleResDOFs.add(resDOFs);
        }
        
        //now rigid-body strand motions...
        ArrayList<DegreeOfFreedom> strandDOFs = strandMotionDOFs(moveableStrands,flexibleRes);
        confDOFs.addAll(strandDOFs);
        
        //...and perturbations
        //standardize conformations first since we'll record initial resBBState here
        standardizeMutatableRes(allowedAAs, flexibleRes);
        
        ArrayList<Perturbation> perts = dset.makePerturbations(m);//will make pert block here
        confDOFs.addAll(perts);
        
        ArrayList<BBFreeBlock> bfbList = getBBFreeBlocks(freeBBZones,flexibleRes);
        for(BBFreeBlock bfb : bfbList)
            confDOFs.addAll( bfb.getDOFs() );
        
        //OK now make RCs using these DOFs
        for(int pos=0; pos<numPos; pos++){
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            
            ArrayList<DegreeOfFreedom> resDOFs = singleResDOFs.get(pos);
            ArrayList<DegreeOfFreedom> resStrandDOFs = strandDOFsForRes(res,strandDOFs);//and check that all res flex on moving strands!
            
            BBFreeBlock curBFB = getCurBFB(bfbList,res);
            
            boolean wtRotOnly = wtRotOnlyRes.contains(flexibleRes.get(pos));
            if(wtRotOnly && !(wtRots.get(pos)!=null&&allowedAAs.get(pos).size()==1) )
                throw new RuntimeException("ERROR: WT rot only on but residue not single AA type with wild-type rotamer");
            
            PositionConfSpace rcs = new PositionConfSpace(pos, res, resDOFs, allowedAAs.get(pos), contSCFlex,
                    resStrandDOFs, perts, dset.getPertIntervals(), dset.getPertStates(pos), curBFB, useEllipses, 
                    wtRots.get(pos), wtRotOnly);
            posFlex.add(rcs);
                        
            if (useEllipses) {
            	confDOFs.addAll(rcs.getEllipsoidalArray());
            } else {
            	confDOFs.addAll(resDOFs);
            }
            
        }
        
        //DEBUG!!!
        /*PDBFileWriter.writePDBFile(m, "STRUCT1.pdb");
        perts.get(0).apply(5);
        PDBFileWriter.writePDBFile(m, "STRUCT2.pdb");
        perts.get(0).apply(0);
        PDBFileWriter.writePDBFile(m, "STRUCT3.pdb");*/
    }
    
    public ConfSpace(ConfSpace other) {
    	// just make a shallow copy
    	this.m = other.m;
    	this.confDOFs = new ArrayList<>(other.confDOFs);
    	this.mutDOFs = new ArrayList<>(other.mutDOFs);
		this.posFlex = new ArrayList<>(other.posFlex);
		this.flexibleRes = new ArrayList<>(other.flexibleRes);
    	this.numPos = other.numPos;
    	this.useEllipses = other.useEllipses;
    }
    
    private ArrayList<BBFreeBlock> getBBFreeBlocks(ArrayList<String[]> freeBBZones, ArrayList<String> flexibleRes){
        //create a BFB for each (start res, end res) pair.  PDB residue numbers provided.  
        ArrayList<BBFreeBlock> ans = new ArrayList<>();
        
        for(String[] termini : freeBBZones){
            ArrayList<Residue> curBFBRes = resListFromTermini(termini, flexibleRes);
            BBFreeBlock bfb = new BBFreeBlock(curBFBRes);
            ans.add(bfb);
        }
        
        return ans;
    }
    
    private BBFreeBlock getCurBFB(ArrayList<BBFreeBlock> bfbList, Residue res){
        //If res is in one of the BFB's, return that BFB; else return null;
        for(BBFreeBlock bfb : bfbList){
            if(bfb.getResidues().contains(res)){//current res is in this BFB
                return bfb;
            }
        }
        
        return null;//not in a BFB
    }
    
    
    private void makeRingPuckers(ArrayList<ArrayList<String>> allowedAAs, ArrayList<String> flexibleRes){
        //For each residue, create a Pro pucker DOF if Pro is allowed
        
        for(int pos=0; pos<numPos; pos++){
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            
            if(res.template.name.equalsIgnoreCase("PRO"))//the reisdue is a proline
                res.pucker = new ProlinePucker(EnvironmentVars.resTemplates, res);
            else {//see if it can mutate to a proline
                for(String AAType : allowedAAs.get(pos)){
                    if(AAType.equalsIgnoreCase("PRO")){
                        res.pucker = new ProlinePucker(EnvironmentVars.resTemplates, res);
                        break;
                    }
                }
            }
            
            if(res.pucker != null)
                confDOFs.add(res.pucker);
        }
    }
    
    
   private ArrayList<Residue> resListFromTermini(String[] termini, ArrayList<String> flexibleRes){ 
       return m.resListFromTermini(termini, flexibleRes);
   }
   
    
    private ArrayList<DegreeOfFreedom> strandMotionDOFs(ArrayList<String[]> moveableStrands,
            ArrayList<String> flexibleRes){
        //Generate all the strand rotation/translation DOFs,
        //given the termini of the moving strands
        
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        
        for(String[] termini : moveableStrands){
            ArrayList<Residue> curStrandRes = resListFromTermini(termini, flexibleRes);
            MoveableStrand str = new MoveableStrand(curStrandRes);
            ans.addAll(str.getDOFs());
        }
        
        return ans;
    }
    
    
    private ArrayList<DegreeOfFreedom> strandDOFsForRes(Residue res, ArrayList<DegreeOfFreedom> strandDOFs){
        //List the strandDOFs that move res
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        
        for(DegreeOfFreedom dof : strandDOFs){
            MoveableStrand str;
            if(dof instanceof StrandRotation)
                str = ((StrandRotation)dof).getMoveableStrand();
            else//must be translation
                str = ((StrandTranslation)dof).getMoveableStrand();
            
            if(str.getResidues().contains(res))
                ans.add(dof);
        }
        
        return ans;
    }
    
    
    private void standardizeMutatableRes(ArrayList<ArrayList<String>> allowedAAs, ArrayList<String> flexibleRes){
        //"mutate" all mutatable residues to the template version of their residue type
        //this ensures that the non-adjustable DOFs (bond angles, etc.) will be as in the template
        //(for consistency purposes)
        for(int pos=0; pos<numPos; pos++){
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            String resName = res.template.name;
            
            if(EnvironmentVars.resTemplates.getTemplate(resName) != null){
                //mutation to current residue type is possible, i.e., this position is mutatable
                //if(HardCodedResidueInfo.canMutateTo(res.template)){//this messes up for N-term mutations
                
                ResidueTypeDOF mutDOF = mutDOFs.get(pos);
                
                for(String allowedAA : allowedAAs.get(pos)){
                    if(allowedAA.equalsIgnoreCase("PRO")){
                        //mutating to and from PRO alters the sidechain a little, including idealization
                        //(because the sidechain connects in two places, it behaves a little different)
                        //Once we mutate to Pro and back we have a consistent conformation
                        mutDOF.mutateTo("PRO");
                        break;
                    }
                }
                
                // AAO 2016: mutation assumes residue is an amino acid. throws an exception otherwise
                if(HardCodedResidueInfo.hasAminoAcidBB(res) && !res.fullName.startsWith("FOL")) {
                    mutDOF.mutateTo(resName);
                }
            }
        }
    }
    
    static ArrayList<DegreeOfFreedom> mutablePosDOFs(Residue res, ArrayList<String> allowedAAs){
        //mutation and dihedral confDOFs for the specified position
        
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        //assuming for now that this is a single-residue position...if multiple just add the following confDOFs for each residue
        
        //mutation DOF
        ans.add(new ResidueTypeDOF(EnvironmentVars.resTemplates, res));
        
        int maxNumDihedrals = 0;//we need to create enough dihedral confDOFs for the allowed AA type with the most dihedrals
        
        for(String AAType : allowedAAs)
            maxNumDihedrals = Math.max( maxNumDihedrals, EnvironmentVars.resTemplates.numDihedralsForResType(AAType) );
        
        for(int dih=0; dih<maxNumDihedrals; dih++)
            ans.add( new FreeDihedral(res,dih) );
                
        if(res.pucker!=null)
            ans.add(res.pucker);
        
        return ans;
    }
    
    
    public double minimizeEnergy(int[] conf, EnergyFunction efunc, String outputPDBFile){
        //minimize the energy of a conformation, within the DOF bounds indicated by conf (a list of RCs)
        //return the minimized energy
        //if outputPDBFile isn't null, then output the minimized conformation to that file
        
        RCTuple RCs = new RCTuple(conf);
        MoleculeModifierAndScorer energy = new MoleculeModifierAndScorer(efunc,this,RCs);
        
        DoubleMatrix1D optDOFVals;
        
        if(energy.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
            Minimizer min = new CCDMinimizer(energy,false);
            //with the generic objective function interface we can easily include other minimizers though

            
            //DEBUG!!!!!  Timing pre-minimization without PB
            /*ArrayList<EnergyFunction> terms = ((MultiTermEnergyFunction)energy.getEfunc()).getTerms();
            ArrayList<Double> coeffs = ((MultiTermEnergyFunction)energy.getEfunc()).getCoeffs();
            if( terms.get(terms.size()-1) instanceof PoissonBoltzmannEnergy ){
                PoissonBoltzmannEnergy pbe = (PoissonBoltzmannEnergy)terms.remove(terms.size()-1);
                double pbcoeff = coeffs.remove(terms.size()-1);
                long startTime = System.currentTimeMillis();
                DoubleMatrix1D startDOFVals = min.minimize();
                long time1 = System.currentTimeMillis();
                terms.add(pbe);
                coeffs.add(pbcoeff);
                ((CCDMinimizer)min).setInitVals(startDOFVals);
                optDOFVals = min.minimize();
            }
            else//NON-DEBUG!*/
                optDOFVals = min.minimize().dofValues;
        }
        else//molecule is already in the right, rigid conformation
            optDOFVals = DoubleFactory1D.dense.make(0);
        
        double minE = energy.getValue(optDOFVals);//this will put m into the minimized conformation
        
        if (outputPDBFile!=null) {
            PDBIO.writeFile(m, "", minE, outputPDBFile);
        }
        
        return minE;
    }
    
    
	public MultiTermEnergyFunction getDecomposedMinimizedEnergy(int[] conf, EnergyFunction efunc, String outputPDBFile){
		//minimize the energy of a conformation, within the DOF bounds indicated by conf (a list of RCs)
		//return the minimized energy
		//if outputPDBFile isn't null, then output the minimized conformation to that file

		RCTuple RCs = new RCTuple(conf);
		MoleculeModifierAndScorer energy = new MoleculeModifierAndScorer(efunc,this,RCs);

		DoubleMatrix1D optDOFVals;

		if(energy.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
			Minimizer min = new CCDMinimizer(energy,false);
			optDOFVals = min.minimize().dofValues;
		}
		else//molecule is already in the right, rigid conformation
			optDOFVals = DoubleFactory1D.dense.make(0);

		double minE = energy.getValue(optDOFVals);//this will put m into the minimized conformation
		
		if (outputPDBFile!=null) {
			PDBIO.writeFile(m, "", minE, outputPDBFile);
		}
		
		return (MultiTermEnergyFunction)energy.getEfunc();
	}
    
    
    public int[] getNumRCsAtPos(){
        //list number of RCs at each position
        int[] numAllowed = new int[numPos];
        
        for(int pos=0; pos<numPos; pos++)
            numAllowed[pos] = posFlex.get(pos).RCs.size();
        
        return numAllowed;
    }
    
    
    public double getRCResEntropy(int pos, int rc){
        String resType = posFlex.get(pos).RCs.get(rc).AAType;
        double resEntropy = EnvironmentVars.resTemplates.getResEntropy(resType);
        return resEntropy;
    }
    
    public double getConfResEntropy(int[] conf){
        double ans = 0;
        
        for(int pos=0; pos<numPos; pos++)
            ans += getRCResEntropy(pos, conf[pos]);
        
        return ans;
    }


	public BigInteger getNumConformations() {
		BigInteger count = BigInteger.valueOf(1);
		for (int pos=0; pos<numPos; pos++) {
			count = count.multiply(BigInteger.valueOf(posFlex.get(pos).RCs.size()));
		}
		return count;
	}
        
    public ArrayList<DegreeOfFreedom> listAllDOFs(){
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        ans.addAll(mutDOFs);
        ans.addAll(confDOFs);
        return ans;
    }
    
    
    /*DoubleMatrix1D[] convertConfToDOFBounds(int[] conf){
        //Argument: RC assignments for all the flexible residues (RCs defined in resFlex)
        //return bounds (lower bounds, upper bounds) for all the degrees of freedom in the system
        RCTuple fullTuple = new RCTuple(conf);
        
        return RCTupleDOFBounds(fullTuple,confDOFs);
        //PULL FROM RC OBJECTS
    }
    
    public DoubleMatrix1D[] RCTupleDOFBounds(RCTuple tuple, ArrayList<DegreeOfFreedom> dofList){
        //Compute the bounds on the confDOFs in DOFlist imposed by the RCs in tuple
        //Create two vectors, vector 0 giving the lower bound for each DOF in dofList,
        //vector 1 giving the upper bound
        //can leave entries at 0 for AA types
        
        //WE SHOUDL TAKE MUT OUT OF DOFS ITS FUNDAMENTALLY DIFFERENT...NO DOUBLE REP, NOT CONF CHANGE
        
        //THIS IS FOR MAKING MOLEC E OBJ FUNCTION
        //MAKE IT TAKE RCs, RETURN ALL CONTINUOUS DOFS AND THEIR BOUNDS
        
        //ACTUALLY WAIT LETS MAKE THIS MOLEC E OBJ FUNCTION CONSTRUCTOR
    }*/
    
    //pairwise minimization will be similar...just only use degrees of freedom affecting the residue pair
    //and the objective function will represent the pairwise energy between the residues
    
    
    
}
