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

package edu.duke.cs.osprey.dof.deeper;

import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;

import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * Holds parameters that control a DEEPer run
 * 
 * @author mhall44
 */
public class DEEPerSettings implements Serializable {
    
	private static final long serialVersionUID = 1L;

	PertSet perts = null;//Describes the actual perturbations
    
    //various settings from ConfigFileParser
    boolean doPerturbations;//We're using DEEPer
    String pertFileName;//The file describing our perturbations

    
    //things for selecting perturbations
    boolean selectPerturbations;//this flag is only used if doPerturbations is true
    String startingPertFile;
    boolean onlyStarting;
    double maxShearParam;
    double maxBackrubParam;
    boolean selectLCAs;
    boolean doRamaCheck;
    //We will need the following basic information about the design system
    ArrayList<String> flexibleRes;
    String PDBFile;
    
    ResidueTemplateLibrary templateLib;

    
    public DEEPerSettings(){
        //Default: no DEEPer
        doPerturbations = false;
    }
    
    public DEEPerSettings(boolean doPerturbations, String pertFileName,
            boolean selectPerturbations, String startingPertFile, boolean onlyStarting, 
            double maxShearParam, double maxBackrubParam, boolean selectLCAs, 
            ArrayList<String> flexibleRes, String PDBFile, boolean doRamaCheck,
						  ResidueTemplateLibrary templateLib) {
        
        this.doPerturbations = doPerturbations;
        this.pertFileName = pertFileName;
        this.selectPerturbations = selectPerturbations;
        this.startingPertFile = startingPertFile;
        this.onlyStarting = onlyStarting;
        this.maxShearParam = maxShearParam;
        this.maxBackrubParam = maxBackrubParam;
        this.selectLCAs = selectLCAs;
        this.flexibleRes = flexibleRes;
        this.PDBFile = PDBFile;
        this.doRamaCheck = doRamaCheck;
        this.templateLib = templateLib;
    }
    
    public boolean doPerturbations() {
    	return doPerturbations;
    }
    
    public void loadPertFile(ResidueTermini termini){
        //load the perturbation file; select perturbations if there is none
        
        if(!doPerturbations)//No perturbations: leave perts null
            return;
        
        perts = new PertSet();
        
        if( !perts.loadPertFile(pertFileName, true, termini) ){//file not found
            
            if(!selectPerturbations)
                throw new RuntimeException("ERROR: Perturbation file not found"
                        + " but not supposed to select perturbations");
            
            PerturbationSelector sele = new PerturbationSelector(startingPertFile, onlyStarting, 
                    maxShearParam, maxBackrubParam, selectLCAs, flexibleRes, PDBFile, termini, 
                    doRamaCheck, templateLib);
            
            PertSet ps = sele.selectPerturbations(termini);
            ps.writePertFile(pertFileName);
            perts.loadPertFile(pertFileName, true, termini);
        }
    }
    
    
    public ArrayList<Perturbation> makePerturbations(Molecule m){
        //Implement the perturbations described in perts in a molecule
        if(perts==null)//no perturbations
            return new ArrayList<>();
        else
            return perts.makePerturbations(m);
    }
    
    
    public ArrayList<ArrayList<double[]>> getPertIntervals(){
        if(perts==null)
            return null;
        
        return perts.pertIntervals;
    }
    
    public ArrayList<ArrayList<int[]>> getPertStates(int pos){
        if(perts==null || perts.pertStates.isEmpty())
            return null;
        
        return perts.pertStates.get(pos);
    }
    
    
    public DEEPerSettings makeDiscreteVersion(){
        //make a version of these settings with continuous intervals changed to discrete points
        DEEPerSettings discrSettings = new DEEPerSettings(doPerturbations, pertFileName+".DISCR", 
            selectPerturbations, startingPertFile, onlyStarting, 0, 0, 
            selectLCAs, flexibleRes, PDBFile, doRamaCheck, templateLib);
        
        if(perts==null)
            discrSettings.perts = null;
        else
            discrSettings.perts = perts.makeDiscreteVersion();
        
        return discrSettings;
    }
    
}
