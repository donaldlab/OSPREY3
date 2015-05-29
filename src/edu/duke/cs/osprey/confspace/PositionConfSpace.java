/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.EllipseCoordDOF;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Residue;

import java.util.ArrayList;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;

/**
 *
 * @author mhall44
 */
public class PositionConfSpace {
    //This class defines the conformational space of a flexible residue
    //including allowed amino-acid types, and rotamers/RCs for each type 
    //subclass PositionConfSpace to make super-residues with super-RCs...
    
    
    public ArrayList<RC> RCs = new ArrayList<>();
    
    public Residue res;//The residue involved
    
    static double dihedFlexInterval = 9;// +/- 9 degree sidechain dihedral continuous flexibility...
    //later can allow this to vary across different dihedrals
    
    ArrayList<EllipseCoordDOF> ellipsoidalDOFs;
    
    public PositionConfSpace(
    		Residue res,
    		ArrayList<DegreeOfFreedom> resDOFs,
    		ArrayList<String> allowedAAs,
    		boolean contSCFlex,
    		boolean useEllipses) {
    	this.res = res;
    	if (useEllipses) { addEllipsoidalRCs(resDOFs, allowedAAs, contSCFlex); } 
    	else { addDihedralRCs(resDOFs, allowedAAs, contSCFlex); }
    }
    
    private void addDihedralRCs(ArrayList<DegreeOfFreedom> resDOFs, ArrayList<String> allowedAAs, boolean contSCFlex) {
        ResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
        
        for(String AAType : allowedAAs){
        	System.out.println("Loading AAType "+AAType);
            int numDihedrals = templateLib.numDihedralsForResType(AAType);
            int numRot = templateLib.numRotForResType(AAType);
                        
            //resDOFs is all sidechain DOFs, for now
            ArrayList<DegreeOfFreedom> dofListForRC = new ArrayList<>();
            for(int dih=0; dih<numDihedrals; dih++) {//get the first numDihedrals dihedrals
                dofListForRC.add(resDOFs.get(dih));
                FreeDihedral x = (FreeDihedral) resDOFs.get(dih);
                System.out.println("Residue: "+x.getResidue().fullName);
            }
                        
            if(numRot==0){//ALA or GLY: no rotamers or dihedrals, so create a single rigid RC
                RC newRC = new RC(AAType, -1, dofListForRC, new ArrayList<Double>(), new ArrayList<Double>(), RCs.size());
                RCs.add(newRC);
            }
                
            for(int rot=0; rot<numRot; rot++){
                //create RC
                ArrayList<Double> dofLB = new ArrayList<>();//lower bounds on each DOF for this RC
                ArrayList<Double> dofUB = new ArrayList<>();//upper bounds

                for(int dih=0; dih<numDihedrals; dih++){
                    double dihedralValForRot = templateLib.getDihedralForRotamer(AAType,rot,dih);
                    
                    if(contSCFlex){//allow continuous flexibility up to dihedFlexInterval in each direction
                        dofLB.add(dihedralValForRot-dihedFlexInterval);
                        dofUB.add(dihedralValForRot+dihedFlexInterval);
                    }
                    else {
                        dofLB.add(dihedralValForRot);
                        dofUB.add(dihedralValForRot);
                    }
                }
                
                RC newRC = new RC(AAType, rot, dofListForRC, dofLB, dofUB, RCs.size());
                RCs.add(newRC);
            }
        }
    	
    }
       
    private void addEllipsoidalRCs(ArrayList<DegreeOfFreedom> resDOFs, ArrayList<String> allowedAAs, boolean contSCFlex) {
    	// resDOFs is just a list of n dihedral angles
    	ResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
    	
    	for (String AAType : allowedAAs) {
    		int numDihedrals = templateLib.numDihedralsForResType(AAType);
    		int numRot = templateLib.numRotForResType(AAType);

    		ArrayList<DegreeOfFreedom> dofListForRC = new ArrayList<>();

    		for (int dih=0; dih<numDihedrals; dih++) { dofListForRC.add(resDOFs.get(dih)); }
    		
    		if (numRot==0) { //if GLY or ALA, just put in a standard rigid rotamer - don't bother with ellipsing
            	ellipsoidalDOFs = new ArrayList<EllipseCoordDOF>();
    			RC newRC = new RC(AAType, -1, dofListForRC, new ArrayList<Double>(), new ArrayList<Double>(), RCs.size());
    			RCs.add(newRC);
    		}
    		
    		for (int rot=0; rot<numRot; rot++) {
    			ArrayList<Double> dofLB = new ArrayList<>();
    			ArrayList<Double> dofUB = new ArrayList<>();
    			
    			double[] dihValues = new double[numDihedrals];
    			for (int dih=0; dih<numDihedrals; dih++) {
    				dihValues[dih] = templateLib.getDihedralForRotamer(AAType, rot, dih);
    			}
    			
    			// generate the list of ellipsoidal DOFs
    			double[] ellValues = getEllipsoidalCoords(dihValues);
    			ArrayList<EllipseCoordDOF> ellCoords = new ArrayList<EllipseCoordDOF>();
    	    	DoubleMatrix2D A = DoubleFactory2D.dense.identity(ellValues.length);
    	    	DoubleMatrix1D c = DoubleFactory1D.dense.make(new double[ellValues.length]);
    	    	for (int i=0; i<ellValues.length; i++) {
    	    		EllipseCoordDOF ellDOF = new EllipseCoordDOF(
    	    				(i==0),
    	    				i,
    	    				ellValues[i],
    	    				A,
    	    				c,
    	    				dofListForRC,
    	    				dihValues);
    	    		ellCoords.add(ellDOF);
    	    		
    	    	}
    	    	
    	    	for (int i=0; i<numDihedrals; i++) {
    	    		if (contSCFlex) {
    	    			dofLB.add(0.0);
    	    			dofUB.add((i==0) ? 300 : 
    	    					(i==numDihedrals-1) ? 2*Math.PI : Math.PI);
    	    		} else {
    	    			dofLB.add(ellCoords.get(i).getCurVal());
    	    			dofUB.add(ellCoords.get(i).getCurVal());
    	    		}
    	    	}
    	    	
    	    	ArrayList<DegreeOfFreedom> dofList = new ArrayList<>();
    	    	dofList.addAll(ellCoords);
    	    	ellipsoidalDOFs = ellCoords;
    	    	RC newRC = new RC(AAType, rot, dofList, dofLB, dofUB, RCs.size());
    	    	RCs.add(newRC);    			
    		}    		
    	}
    }
    
	public double[] getEllipsoidalCoords(double[] dihedrals) {
    	if (dihedrals.length==0) { return new double[0]; }
    	
    	// for now we're just using the unit sphere
    	DoubleMatrix2D A = DoubleFactory2D.dense.identity(dihedrals.length);
    	DoubleMatrix1D c = DoubleFactory1D.dense.make(new double[dihedrals.length]);
    	
		EigenvalueDecomposition evd = new EigenvalueDecomposition(A);
		DoubleMatrix2D Q = evd.getV();
		DoubleMatrix2D L = evd.getD();
		DoubleMatrix2D qT = Q.viewDice().copy();
		Algebra alg = new Algebra();
    	
    	// first transform the cartesian coordinates based on the ellipse
    	double[] s = new double[dihedrals.length];
    	for (int i=0; i<dihedrals.length; i++) { 
    		s[i] = dihedrals[i]-c.get(i);
    	}
    	double[] u = alg.mult(qT, DoubleFactory1D.dense.make(s)).toArray();
    	double[] x = new double[u.length];
    	for (int i=0; i<u.length; i++) {
    		x[i] = u[i]/Math.sqrt(L.get(i, i));
    	}    	
    	dihedrals = x;
    	
    	// now get elliptical coordinates
    	double radius = 0;
    	for (double d : dihedrals) { radius += d*d; }
    	radius = Math.sqrt(radius);
    	int n = dihedrals.length;
    	double[] phi = new double[n-1];
    	for (int i=0; i<n-1; i++) {
    		double d=0;
    		for (int j=i; j<n; j++) { d += dihedrals[j]*dihedrals[j]; }
    		double quot = dihedrals[i]/Math.sqrt(d);
    		phi[i] = Math.acos(quot);
    	}
    	if (dihedrals[n-1] < 0) { phi[n-2] = 2*Math.PI - phi[n-2]; }
    	double[] ellCoords = new double[n];
    	ellCoords[0] = radius;
    	for (int i=1; i<n; i++) { ellCoords[i] = phi[i-1]; }
    	return ellCoords;
    }

    public ArrayList<EllipseCoordDOF> getEllipsoidalArray() {
    	return this.ellipsoidalDOFs;
    }
}
























