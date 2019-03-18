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

package edu.duke.cs.osprey.dof;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.TreeSet;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RigidBodyMotion;

/**
 *
 * @author mhall44
 */

// TODO: this isn't really a DOF
public class ResidueTypeDOF extends DegreeOfFreedom {
    //This degree of freedom is the residue type (e.g., AA type) at a particular position
    //So applying values of it means mutating the residue
    
    private static final long serialVersionUID = 4285811771185813789L;
    
    // TODO: this should be final and not transient, but we're stuck using weird serialization for now
    public transient ResidueTemplateLibrary templateLib;
    
    private Residue res;//what residue in the molecule we are talking about
    private boolean idealizeSidechainAfterMutation;
    
    /**
     * this doesn't need to be a DOF anymore, just the static function switchToTemplate()
     */
    @Deprecated
    public ResidueTypeDOF(ResidueTemplateLibrary templateLib, Residue res) {
        this(templateLib, res, false);
    }
    
    /**
     * this doesn't need to be a DOF anymore, just the static function switchToTemplate()
     */
    @Deprecated
    public ResidueTypeDOF(ResidueTemplateLibrary templateLib, Residue res, boolean idealizeSidechainAfterMutation) {
        this.templateLib = templateLib;
        this.res = res;
        this.idealizeSidechainAfterMutation = idealizeSidechainAfterMutation;
    }
    
    // TEMP TODO HACKHACK: Java's serialization system forces us to de-serialize object without any context
    // but these DoFs need the template library, which apparently doesn't serialize correctly
    // (it's stupid... we shouldn't need to serialize the template library as part of the EPIC matrix anyway...)
    // so explicitly force Java's default de-serializer to use the EnvironmentVars
    // for now... need to find a better way to do this in the future
    private void readObject(ObjectInputStream ois)
    throws ClassNotFoundException, IOException {
        ois.defaultReadObject();
        templateLib = EnvironmentVars.resTemplates;
        assert (templateLib != null);
    }
    
    public void mutateTo(String resType) {
        switchToTemplate(getLibraryTemplate(resType));
    }
    
    public ResidueTemplate getLibraryTemplate(String resType) {
        return templateLib.getTemplateForMutation(resType, res);
    }
    
    public boolean isTemplate(ResidueTemplate template) {
    	return this.res.template == template;
    }
    
    public void switchToTemplate(ResidueTemplate newTemplate) {
    	switchToTemplate(templateLib, res, newTemplate, idealizeSidechainAfterMutation);
    }

	public static void switchToTemplate(ResidueTemplateLibrary templateLib, Residue res, ResidueTemplate newTemplate, boolean idealizeSidechainAfterMutation) {
    	switchToTemplate(templateLib, res, newTemplate, idealizeSidechainAfterMutation, new MutAlignmentCache(), true);
	}

    public static void switchToTemplate(ResidueTemplateLibrary templateLib, Residue res, ResidueTemplate newTemplate, boolean idealizeSidechainAfterMutation, MutAlignmentCache mutAlignmentCache, boolean reconnectInterResBonds) {
        ResidueTemplate oldTemplate = res.template;

        if (oldTemplate.CAEquivalent == null || newTemplate.CAEquivalent == null) {//non-mutatable templates
            if (oldTemplate.name.equalsIgnoreCase(newTemplate.name))//let it be so we can make non-mutatable residue types flexible
                return;
            else {
                throw new RuntimeException("ERROR: Trying to mutate " + oldTemplate.name + " to " + newTemplate.name +
                        " but CAEQUIVALENT not specified in template and cannot be inferred");
            }
        }

        //the residue's going to change some, so break its inter-residue bonds
        res.removeInterResBonds();
        res.intraResBondsMarked = false;//we'll need to redo these too
        
        res.template = newTemplate;
        
        res.fullName = newTemplate.name + res.fullName.substring(3);
        //res type name is first three characters of full name
        
        MutAlignment mutAlignment = mutAlignmentCache.get(oldTemplate, newTemplate);
        //coordinates will come from the template,
        //but we'll move them as a rigid body to match the backbone atoms
        int[][] mutAlignAtoms = mutAlignment.getMutAlignmentAtoms();
        //2x4 array of which atoms are used in the old and new residues to align
        //to do the mutation
        
        double oldMAACoords[][] = extractCoords(mutAlignAtoms[0],res.coords);
        //coordinates of the atoms in the old residue to align
        
        double newCoords[] = newTemplate.templateRes.coords.clone();
        double templateMAACords[][] = extractCoords(mutAlignAtoms[1],newCoords);
        
        //we now construct a rigid-body motion that will map the sidechain (or generally,
        //the non-"backbone" part, if not a standard amino acid) from the template to the residue's
        //curent frame of reference
        RigidBodyMotion templMotion = new RigidBodyMotion(templateMAACords,oldMAACoords);
        
        templMotion.transform(newCoords);
        
        //the backbone atoms will be kept exactly as before the mutation
        //if the sidechain attaches only to the first mutAlignAtom, this method keeps bond lengths
        //exactly as in the template for sidechain, and as in the old backbone otherwise
        TreeSet<String> BBAtomNames =  mutAlignment.getNonmovingAtomNames();
        for(String BBAtomName : BBAtomNames){
            int BBAtomIndexOld = oldTemplate.templateRes.getAtomIndexByName(BBAtomName);
            int BBAtomIndexNew = newTemplate.templateRes.getAtomIndexByName(BBAtomName);
            
            //copy coordinates of the BB atom from old to new coordinates
            if(BBAtomIndexOld!=-1 && BBAtomIndexNew!=-1)
                System.arraycopy(res.coords, 3*BBAtomIndexOld, newCoords, 3*BBAtomIndexNew, 3);
        }
        
        res.coords = newCoords;
        
        //finally, update atoms in res to match new template
        ArrayList<Atom> newAtoms = new ArrayList<>();
        for(Atom at : newTemplate.templateRes.atoms){
            Atom newAtom = at.copy();
            newAtom.res = res;
            newAtoms.add(newAtom);
        }
        res.atoms = newAtoms;    
        
        //reconnect all bonds
        res.markIntraResBondsByTemplate();
        if (reconnectInterResBonds) {
        	res.reconnectInterResBonds();
		}
        
        //special case if sidechain loops back in additional place to backbone...
        if(oldTemplate.name.equalsIgnoreCase("PRO") || newTemplate.name.equalsIgnoreCase("PRO")){
            if (idealizeSidechainAfterMutation) {
                SidechainIdealizer.idealizeSidechain(templateLib, res);
            }
            if(!newTemplate.name.equalsIgnoreCase("PRO")){//if mutating from Pro, no ring closure issues possible anymore
                if(res.pucker!=null){
                    if(res.pucker.puckerProblem != null){
                        res.pucker.puckerProblem.removeFromRes();
                        res.pucker.puckerProblem = null;
                    }
                }
            }
        }
        else if(idealizeSidechainAfterMutation && HardCodedResidueInfo.hasAminoAcidBB(res)
                && (res.template.name.equalsIgnoreCase("GLY")||res.getAtomByName("CB")!=null)){//trying to make sure it's really an amino acid
            SidechainIdealizer.idealizeSidechain(templateLib, res);
        }
    }
    
    public void restoreCoordsFromTemplate() {

        if(res.template.CAEquivalent==null)//can't do this for non-mutable residues
            return;
    
        // get the alignment of backbone atoms
        MutAlignment mutAlignment = new MutAlignment(res.template, res.template);
        int[][] mutAlignAtoms = mutAlignment.getMutAlignmentAtoms();
        double resBBCoords[][] = extractCoords(mutAlignAtoms[0], res.coords);
        double templateBBCoords[][] = extractCoords(mutAlignAtoms[1], res.template.templateRes.coords);
        
        // rotation from template to res
        RigidBodyMotion xform = new RigidBodyMotion(templateBBCoords, resBBCoords);
        TreeSet<String> nonmovingAtomNames = mutAlignment.getNonmovingAtomNames();
        for (Atom atom : res.atoms) {
            
            // skip backbone atoms
            if (nonmovingAtomNames.contains(atom.name)) {
                continue;
            }
            
            // transform sidechain atoms
            int i = atom.indexInRes*3;
            System.arraycopy(res.template.templateRes.coords, i, res.coords, i, 3);
            xform.transform(res.coords, atom.indexInRes);
        }
    }
    
    
    static double[][] extractCoords(int[] index, double[] allCoords){
        //allCoords is the concatenated coords of a bunch of atoms
        //get the coords (3-D) of the atoms with the specified indices
        double[][] ans = new double[index.length][3];
        
        for(int i=0; i<index.length; i++){
            System.arraycopy(allCoords, 3*index[i], ans[i], 0, 3);
        }
        
        return ans;
    }
    
    
    public String getCurResType(){//current residue type
        return res.fullName.substring(0,3);//this is always the first three letters of the full name
    }
    

    @Override
    public void apply(double paramVal) {
        throw new IllegalArgumentException("ERROR: ResidueTypeDOF takes an residue type name"
                + " as argument; can't take "+paramVal);
    }
    
    
    
    @Override
    public Residue getResidue() { return res; }
    
    
    @Override
    public DegreeOfFreedom copy() {
        return new ResidueTypeDOF(templateLib, res);
    }
    
    @Override
    public void setMolecule(Molecule val) {
        res = val.getResByPDBResNumber(res.getPDBResNumber());
    }

    @Override
    public DOFBlock getBlock(){
        return null;
    }

    @Override
    public String getName() {
        return "RESTYPE"+res.getPDBResNumber();
    }
}
