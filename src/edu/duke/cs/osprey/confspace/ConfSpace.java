/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.EllipseCoordDOF;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MolecEObjFunction;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;

import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class ConfSpace {
    //This class represents the conformational search space for a design
    //used for GMEC-based design, K*, or anything else we want to do with a conformational space
    //This class just defines the conformational space itself (the molecule + all kinds of flexibility
    //and possible mutations, and how these are represented as RCs, etc.)
    //This class can be put in an AnnotatedConfSpace to add annotations like what RCs are pruned,
    //what their pairwise energies are, etc.  
    
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

    
    public int numPos;//number of flexible positions
    
    public boolean useEllipses = false;
    
    
    public ConfSpace(
    		String PDBFile,
    		ArrayList<String> flexibleRes,
    		ArrayList<ArrayList<String>> allowedAAs,
    		boolean addWT,
            boolean contSCFlex,
            boolean ellipses){
        //initialize a new conformational space, desomefining all its flexibility
        //we use one residue per position here
        //PDBFile: the structure to read from
        //flexibleRes: list of residue numbers to be made flexible (as in PDB file)
        //allowedAAs: list of allowed residue types at each flexible position
        //addWT: whether to add wild-type to the allowed AA types
        
        //FLEXIBILITY: We do a rotamer search in all cases
        //contSCFlex means allow continuous sidechain flexibility
        //ADD OTHER OPTIONS: WT ROTAMERS, DIFFERENT ROT WIDTHS, DEEPER, RIGID-BODY MOTIONS
        
    	useEllipses = ellipses;  	
    	
        numPos = flexibleRes.size();
        
        //read the structure and assign templates, deleting unassignable res...
        m = PDBFileReader.readPDBFile(PDBFile);
                
        for(int pos=0; pos<numPos; pos++){
            
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            
            if(addWT){//at this point, m has all wild-type residues, so just see what res is now
                String wtName = res.template.name;
                if(!allowedAAs.get(pos).contains(wtName))
                    allowedAAs.get(pos).add(wtName);
            }
            
            ArrayList<DegreeOfFreedom> resDOFs = mutablePosDOFs(res,allowedAAs.get(pos));//add mutation and dihedral confDOFs
                        
            ResidueTypeDOF resMutDOF = (ResidueTypeDOF)resDOFs.remove(0);//first mutable pos DOF is the mutation-type DOF
            
            if (useEllipses) {
            	if (resDOFs.size() == 0) { // empty ArrayList if no DOFs
            		confDOFs.addAll(new ArrayList<EllipseCoordDOF>());
            	} else {   	
	            	// generate radius
	            	double radius = 0;
	            	for (DegreeOfFreedom resDOF : resDOFs) { 
	            		radius += resDOF.getCurVal()*resDOF.getCurVal();
	            	}
	            	radius = Math.sqrt(radius);            	
	            	
	            	// generate angles
	            	double[] angles = new double[resDOFs.size()-1];
	            	int n = angles.length;
	            	for (int i=0; i<n; i++) {
	            		double diff = 0;
	            		for (int j=0; j<i; j++) { diff += Math.pow(resDOFs.get(j).getCurVal(), 2); }
	            		double den = Math.sqrt(radius*radius - diff);
	            		angles[i]=Math.acos(resDOFs.get(i).getCurVal()/den);
	            	}            	
	            	            	
	            	// for now, just using a basic sphere with radius 0.5 radians
	            	double[] center = new double[n+1];
	            	for (int i=0; i<resDOFs.size(); i++) { 
	            		center[i] = resDOFs.get(i).getCurVal();
	            	}
	            	DoubleMatrix1D c = DoubleFactory1D.dense.make(center);
	            	DoubleMatrix2D a = DoubleFactory2D.dense.identity(n+1);
	            	// TODO:
	            	// 		identify what rotamer each RC is at to get the right ellipse
	            	// 		load the matrices from some other data construct (i.e. a class)
	            	// 		perhaps use the dihedrals to identify a point, and derive the
	            	// 		cluster from what surrounds that point? 
	            	
	            	// make the EllipseCoordDOF array and add it to ConfDOFs
	            	ArrayList<EllipseCoordDOF> ellipseDOFs = new ArrayList<EllipseCoordDOF>();
	            	ellipseDOFs.add(new EllipseCoordDOF(true, 0, a, c, resDOFs)); //radius
	            	for (int i=0; i<angles.length; i++) {
	            		ellipseDOFs.add(new EllipseCoordDOF(false, angles[i], a, c, resDOFs));
	            	}
	            	confDOFs.addAll(ellipseDOFs);
            	}
            } else { 
            	confDOFs.addAll(resDOFs);//record the conformational confDOFs here
            }
            mutDOFs.add(resMutDOF);
            
            PositionConfSpace rcs = new PositionConfSpace(res,resDOFs,allowedAAs.get(pos),contSCFlex);
            posFlex.add(rcs);
        }
        
        standardizeMutatableRes();
    }
    
    
    private void standardizeMutatableRes(){
        //"mutate" all mutatable residues to the template version of their residue type
        //this ensures that the non-adjustable DOFs (bond angles, etc.) will be as in the template
        //(for consistency purposes)
        for(int pos=0; pos<numPos; pos++){
            Residue res = posFlex.get(pos).res;
            if(HardCodedResidueInfo.canMutateTo(res.template)){
                ResidueTypeDOF mutDOF = mutDOFs.get(pos);
                mutDOF.mutateTo(res.template.name);
            }
        }
    }
    
    static ArrayList<DegreeOfFreedom> mutablePosDOFs(Residue res, ArrayList<String> allowedAAs){//mutation and dihedral confDOFs for the specified position
        
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        //assuming for now that this is a single-residue position...if multiple just add the following confDOFs for each residue
        
        //mutation DOF
        ans.add(new ResidueTypeDOF(res));
        
        int maxNumDihedrals = 0;//we need to create enough dihedral confDOFs for the allowed AA type with the most dihedrals
        
        for(String AAType : allowedAAs)
            maxNumDihedrals = Math.max( maxNumDihedrals, EnvironmentVars.resTemplates.numDihedralsForResType(AAType) );
        
        for(int dih=0; dih<maxNumDihedrals; dih++)
            ans.add( new FreeDihedral(res,dih) );
        
        return ans;
    }
    
    
    public double minimizeEnergy(int[] conf, EnergyFunction efunc, String outputPDBFile){
        //minimize the energy of a conformation, within the DOF bounds indicated by conf (a list of RCs)
        //return the minimized energy
        //if outputPDBFile isn't null, then output the minimized conformation to that file
        
        RCTuple RCs = new RCTuple(conf);
        MolecEObjFunction energy = new MolecEObjFunction(efunc,this,RCs);
        
        DoubleMatrix1D optDOFVals;
        
        if(energy.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
            Minimizer min = new CCDMinimizer(energy,false);
            //with the generic objective function interface we can easily include other minimizers though

            optDOFVals = min.minimize();
        }
        else//molecule is already in the right, rigid conformation
            optDOFVals = DoubleFactory1D.dense.make(0);
        
        double minE = energy.getValue(optDOFVals);//this will put m into the minimized conformation
        
        if(outputPDBFile!=null)
            PDBFileWriter.writePDBFile(m, outputPDBFile, minE);
        
        return minE;
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
